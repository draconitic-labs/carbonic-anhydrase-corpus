#!/usr/bin/env python3
"""Build a deterministic human alpha-carbonic-anhydrase evidence corpus.

Only official UniProt, Rhea, PDBe/SIFTS, and ChEMBL APIs are used. Network
responses are retained byte-for-byte under ``raw/`` with retrieval metadata and
checksums. A normal run is cache-first; use ``--refresh`` to replace snapshots
or ``--offline`` to prove the normalized corpus can be reproduced from cache.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = "1.0.0"
FAMILY_ID = "human-alpha-carbonic-anhydrase"
USER_AGENT = "Draconitic-evidence-corpus/1.0 (+https://draconitic.com)"

# Accession pins are independently checked against each live/cached UniProt
# entry. Quotas sum to the declared release counts (80 / 28 / 12).
CORE_PROTEINS: tuple[dict[str, Any], ...] = (
    {"gene": "CA1", "accession": "P00915", "catalytic": True, "structures": 10, "activities": 2},
    {"gene": "CA2", "accession": "P00918", "catalytic": True, "structures": 20, "activities": 4},
    {"gene": "CA3", "accession": "P07451", "catalytic": True, "structures": 5, "activities": 2},
    {"gene": "CA4", "accession": "P22748", "catalytic": True, "structures": 10, "activities": 2},
    {"gene": "CA5A", "accession": "P35218", "catalytic": True, "structures": 0, "activities": 2},
    {"gene": "CA5B", "accession": "Q9Y2D0", "catalytic": True, "structures": 0, "activities": 2},
    {"gene": "CA6", "accession": "P23280", "catalytic": True, "structures": 1, "activities": 2},
    {"gene": "CA7", "accession": "P43166", "catalytic": True, "structures": 8, "activities": 2},
    {"gene": "CA8", "accession": "P35219", "catalytic": False, "structures": 1, "activities": 0},
    {"gene": "CA9", "accession": "Q16790", "catalytic": True, "structures": 8, "activities": 2},
    {"gene": "CA10", "accession": "Q9NS85", "catalytic": False, "structures": 0, "activities": 0},
    {"gene": "CA11", "accession": "O75493", "catalytic": False, "structures": 0, "activities": 2},
    {"gene": "CA12", "accession": "O43570", "catalytic": True, "structures": 8, "activities": 2},
    {"gene": "CA13", "accession": "Q8N1Q1", "catalytic": True, "structures": 7, "activities": 2},
    {"gene": "CA14", "accession": "Q9ULX7", "catalytic": True, "structures": 2, "activities": 2},
)

EXPECTED_STRUCTURE_TOTAL = 80
EXPECTED_ACTIVITY_TOTAL = 28
EXPECTED_REACTION_TOTAL = 12
RHEA_MASTER_ID = "RHEA:10748"

LICENSE_SOURCES = (
    (
        "raw/licenses/uniprot-license.txt",
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/LICENSE",
        "text/plain",
    ),
    (
        "raw/licenses/rhea-license.txt",
        "https://ftp.expasy.org/databases/rhea/LICENSE.txt",
        "text/plain",
    ),
)

# Web pages that state a source's terms are kept as a short text extract (the sentence that
# names the licence, the URL, the retrieval time and the page's SHA-256), not as the whole page.
LICENSE_PAGE_EXTRACTS = (
    (
        "raw/licenses/chembl-license-extract.txt",
        "https://chembl.gitbook.io/chembl-interface-documentation/about",
        "Creative Commons Attribution-Share Alike 3.0",
    ),
    (
        "raw/licenses/pdbe-license-extract.txt",
        "https://www.ebi.ac.uk/pdbe/about/public-data-access-statement",
        "Creative Commons Zero",
    ),
)

LICENSE_URLS = {
    "CC-BY-4.0": "https://creativecommons.org/licenses/by/4.0/legalcode",
    "CC-BY-SA-3.0": "https://creativecommons.org/licenses/by-sa/3.0/legalcode",
    "CC0-1.0": "https://creativecommons.org/publicdomain/zero/1.0/legalcode",
}

# The licence that governs each record type: the most restrictive of the sources it copies from.
# Structure mappings copy UniProt cross-references (CC BY 4.0) beside CC0 PDBe data; activities copy
# ChEMBL records (CC BY-SA 3.0); reaction assignments copy UniProt and Rhea (CC BY 4.0).
RECORD_LICENSES = {
    "experimental_structure_mapping": "CC-BY-4.0",
    "assay_activity": "CC-BY-SA-3.0",
    "reaction_assignment": "CC-BY-4.0",
}

def chembl_citation(chembl_release: str) -> str:
    """ChEMBL's requested citation for a release identifier such as ChEMBL_37."""
    number = chembl_release.rsplit("_", 1)[-1]
    return f"ChEMBL release {number}, doi:10.6019/CHEMBL.database.{number}"


def license_page_extract(page: bytes, url: str, retrieved_at: str, phrase: str) -> str:
    """Return the text extract for a licence page: the sentence containing ``phrase``."""
    text = page.decode("utf-8", "replace")
    text = re.sub(r"(?is)<(script|style)\b.*?</\1>", " ", text)
    text = html.unescape(re.sub(r"(?s)<[^>]+>", " ", text))
    text = re.sub(r"\s+", " ", text)
    at = text.find(phrase)
    if at == -1:
        raise RuntimeError(f"Licence statement not found on {url}")
    start = text.rfind(". ", 0, at) + 2 if text.rfind(". ", 0, at) != -1 else 0
    end = text.find(". ", at)
    sentence = text[start : end + 1 if end != -1 else len(text)].strip()
    return (
        f"Extract of {url}\n"
        f"Retrieved {retrieved_at}\n"
        f"Page SHA-256 {sha256_bytes(page)}\n"
        f"\n{sentence}\n"
    )


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    path.write_text(payload, encoding="utf-8")


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")


class SnapshotDownloader:
    """Cache-first downloader that records exact response provenance."""

    def __init__(self, *, refresh: bool, offline: bool) -> None:
        self.refresh = refresh
        self.offline = offline
        self.manifest_path = ROOT / "raw/retrieval-manifest.json"
        self._lock = threading.Lock()
        self._path_locks: dict[str, threading.Lock] = defaultdict(threading.Lock)
        if self.manifest_path.exists():
            loaded = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            self.entries: dict[str, dict[str, Any]] = {
                entry["path"]: entry for entry in loaded.get("entries", [])
            }
        else:
            self.entries = {}

    def fetch(self, rel_path: str, url: str, *, accept: str = "application/json") -> bytes:
        path = ROOT / rel_path
        lock = self._path_locks[rel_path]
        with lock:
            if path.exists() and not self.refresh:
                data = path.read_bytes()
                meta = self.entries.get(rel_path)
                if meta is None:
                    raise RuntimeError(
                        f"Cached snapshot {rel_path} has no retrieval metadata; use --refresh to repair it"
                    )
                if meta["url"] != url:
                    raise RuntimeError(f"Snapshot URL drift for {rel_path}: {meta['url']} != {url}")
                actual = sha256_bytes(data)
                if actual != meta["sha256"]:
                    raise RuntimeError(f"Cached snapshot checksum mismatch: {rel_path}")
                return data

            if self.offline:
                raise RuntimeError(f"Offline build is missing required snapshot: {rel_path}")

            data, headers, status = self._download(url, accept)
            return self._store(rel_path, url, data, headers, status)

    def extract(self, rel_path: str, url: str, phrase: str) -> bytes:
        """Like fetch, but keep only the licence sentence of an HTML page (see license_page_extract)."""
        if (ROOT / rel_path).exists() and not self.refresh:
            return self.fetch(rel_path, url, accept="text/html")
        with self._path_locks[rel_path]:
            if self.offline:
                raise RuntimeError(f"Offline build is missing required snapshot: {rel_path}")
            page, headers, status = self._download(url, "text/html")
            retrieved_at = now_utc()
            data = license_page_extract(page, url, retrieved_at, phrase).encode("utf-8")
            extra = {"page_sha256": sha256_bytes(page), "page_content_length": len(page)}
            return self._store(rel_path, url, data, headers, status, retrieved_at=retrieved_at, extra=extra)

    def _download(self, url: str, accept: str) -> tuple[bytes, dict[str, str], int]:
        request = urllib.request.Request(
            url,
            headers={"Accept": accept, "User-Agent": USER_AGENT},
            method="GET",
        )
        error: Exception | None = None
        for attempt in range(6):
            try:
                with urllib.request.urlopen(request, timeout=90) as response:
                    data = response.read()
                    headers = {key.lower(): value for key, value in response.headers.items()}
                    status = response.status
                return data, headers, status
            except urllib.error.HTTPError as exc:
                error = exc
                if exc.code not in {429, 500, 502, 503, 504}:
                    raise
            except (urllib.error.URLError, TimeoutError) as exc:
                error = exc
            if attempt == 5:
                raise RuntimeError(f"Unable to retrieve {url}: {error}") from error
            time.sleep(min(12, 1.5 * (2**attempt)))
        raise RuntimeError(f"Unable to retrieve {url}: {error}")

    def _store(
        self,
        rel_path: str,
        url: str,
        data: bytes,
        headers: dict[str, str],
        status: int,
        *,
        retrieved_at: str | None = None,
        extra: dict[str, Any] | None = None,
    ) -> bytes:
        path = ROOT / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(data)
        temporary.replace(path)
        selected_headers = {
            key: headers[key]
            for key in ("content-type", "etag", "last-modified", "x-uniprot-release", "x-uniprot-release-date")
            if key in headers
        }
        entry = {
            "path": rel_path,
            "url": url,
            "retrieved_at": retrieved_at or now_utc(),
            "http_status": status,
            "content_length": len(data),
            "sha256": sha256_bytes(data),
            "response_headers": selected_headers,
            **(extra or {}),
        }
        with self._lock:
            self.entries[rel_path] = entry
        return data

    def json(self, rel_path: str, url: str) -> Any:
        return json.loads(self.fetch(rel_path, url, accept="application/json").decode("utf-8"))

    def text(self, rel_path: str, url: str, *, accept: str = "text/plain") -> str:
        return self.fetch(rel_path, url, accept=accept).decode("utf-8")

    def source_ref(self, rel_path: str) -> dict[str, Any]:
        entry = self.entries[rel_path]
        return {
            "snapshot_path": rel_path,
            "snapshot_sha256": entry["sha256"],
            "retrieved_at": entry["retrieved_at"],
            "url": entry["url"],
        }

    def save_manifest(self) -> None:
        write_json(
            self.manifest_path,
            {
                "schema_version": SCHEMA_VERSION,
                "description": "Exact official-source HTTP snapshots used to build the normalized corpus.",
                "entries": [self.entries[path] for path in sorted(self.entries)],
            },
        )


def fetch_json_specs(
    downloader: SnapshotDownloader,
    specs: Iterable[tuple[str, str]],
    *,
    workers: int = 8,
) -> dict[str, Any]:
    unique = {rel_path: url for rel_path, url in specs}
    results: dict[str, Any] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(downloader.json, rel_path, url): rel_path
            for rel_path, url in unique.items()
        }
        for future in as_completed(futures):
            rel_path = futures[future]
            results[rel_path] = future.result()
    return results


def primary_gene(entry: dict[str, Any]) -> str | None:
    for gene in entry.get("genes", []):
        name = gene.get("geneName", {}).get("value")
        if name:
            return name
    return None


def protein_name(entry: dict[str, Any]) -> str | None:
    description = entry.get("proteinDescription", {})
    recommended = description.get("recommendedName", {}).get("fullName", {}).get("value")
    if recommended:
        return recommended
    submissions = description.get("submissionNames") or []
    return submissions[0].get("fullName", {}).get("value") if submissions else None


def xrefs(entry: dict[str, Any], database: str) -> list[dict[str, Any]]:
    return [item for item in entry.get("uniProtKBCrossReferences", []) if item.get("database") == database]


def comment_texts(entry: dict[str, Any], comment_type: str) -> list[str]:
    values: list[str] = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") != comment_type:
            continue
        for text in comment.get("texts", []):
            value = text.get("value")
            if value:
                values.append(value)
    return values


def reaction_with_rhea(entry: dict[str, Any], rhea_id: str) -> dict[str, Any] | None:
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        reaction = comment.get("reaction") or {}
        if any(
            ref.get("database") == "Rhea" and ref.get("id") == rhea_id
            for ref in reaction.get("reactionCrossReferences", [])
        ):
            return reaction
    return None


def property_map(cross_reference: dict[str, Any]) -> dict[str, Any]:
    return {item.get("key"): item.get("value") for item in cross_reference.get("properties", [])}


def parse_resolution(value: Any) -> float | None:
    if not isinstance(value, str):
        return None
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)", value)
    return float(match.group(1)) if match else None


def resolution_band(value: float | None) -> str:
    if value is None:
        return "unreported"
    if value < 1.5:
        return "lt-1.5A"
    if value < 2.0:
        return "1.5-to-lt-2.0A"
    if value < 2.5:
        return "2.0-to-lt-2.5A"
    return "ge-2.5A"


def structure_candidate_order(references: list[dict[str, Any]]) -> list[dict[str, Any]]:
    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for reference in references:
        props = property_map(reference)
        method = str(props.get("Method") or "unreported")
        resolution = parse_resolution(props.get("Resolution"))
        candidate = {
            "reference": reference,
            "method": method,
            "resolution_angstrom": resolution,
            "resolution_band": resolution_band(resolution),
        }
        buckets[f"{method}|{candidate['resolution_band']}"].append(candidate)
    for values in buckets.values():
        values.sort(key=lambda item: item["reference"]["id"])
    queues = {key: deque(values) for key, values in sorted(buckets.items())}
    ordered: list[dict[str, Any]] = []
    while any(queues.values()):
        for key in sorted(queues):
            if queues[key]:
                ordered.append(queues[key].popleft())
    return ordered


def mapping_for_accession(payload: dict[str, Any], pdb_id: str, accession: str) -> dict[str, Any] | None:
    entry = payload.get(pdb_id.lower()) or payload.get(pdb_id.upper()) or {}
    return entry.get("UniProt", {}).get(accession)


def parse_rhea_tsv(payload: str) -> dict[str, str]:
    rows = list(csv.DictReader(io.StringIO(payload), delimiter="\t"))
    matches = [row for row in rows if row.get("Reaction identifier") == RHEA_MASTER_ID]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {RHEA_MASTER_ID} row from Rhea, found {len(matches)}")
    return matches[0]


def normalize_protein(
    cfg: dict[str, Any], entry: dict[str, Any], downloader: SnapshotDownloader
) -> dict[str, Any]:
    gene = cfg["gene"]
    accession = cfg["accession"]
    if entry.get("primaryAccession") != accession:
        raise RuntimeError(f"UniProt accession mismatch for {gene}")
    if primary_gene(entry) != gene:
        raise RuntimeError(f"UniProt primary gene mismatch for {accession}: {primary_gene(entry)}")
    if entry.get("organism", {}).get("taxonId") != 9606:
        raise RuntimeError(f"UniProt taxon mismatch for {accession}")
    if "reviewed" not in str(entry.get("entryType", "")).lower():
        raise RuntimeError(f"UniProt entry is not reviewed: {accession}")

    reaction = reaction_with_rhea(entry, RHEA_MASTER_ID)
    functions = comment_texts(entry, "FUNCTION")
    if cfg["catalytic"] and reaction is None:
        raise RuntimeError(f"Catalytic pin lacks {RHEA_MASTER_ID}: {gene}")
    if not cfg["catalytic"]:
        if reaction is not None:
            raise RuntimeError(f"Acatalytic pin unexpectedly cites {RHEA_MASTER_ID}: {gene}")
        if not any("does not have" in text.lower() and "catalytic" in text.lower() for text in functions):
            raise RuntimeError(f"Acatalytic status is not explicit in UniProt FUNCTION text: {gene}")

    gene_synonyms = []
    for gene_block in entry.get("genes", []):
        gene_synonyms.extend(item.get("value") for item in gene_block.get("synonyms", []) if item.get("value"))

    alternative_names = []
    for name in entry.get("proteinDescription", {}).get("alternativeNames", []):
        value = name.get("fullName", {}).get("value")
        if value:
            alternative_names.append(value)

    pdb_refs = xrefs(entry, "PDB")
    chembl_ids = sorted(item["id"] for item in xrefs(entry, "ChEMBL"))
    sequence = entry.get("sequence", {})
    audit = entry.get("entryAudit", {})
    rel_path = f"raw/uniprot/entries/{accession}.json"
    return {
        "protein_id": f"UNIPROT:{accession}",
        "family_id": FAMILY_ID,
        "gene_symbol": gene,
        "gene_synonyms": sorted(set(gene_synonyms)),
        "uniprot_accession": accession,
        "uniprot_entry_name": entry.get("uniProtkbId"),
        "review_status": "reviewed",
        "organism": entry.get("organism", {}).get("scientificName"),
        "taxon_id": entry.get("organism", {}).get("taxonId"),
        "protein_name": protein_name(entry),
        "alternative_names": alternative_names,
        "sequence": {
            "value": sequence.get("value"),
            "length": sequence.get("length"),
            "molecular_weight": sequence.get("molWeight"),
            "crc64": sequence.get("crc64"),
            "modified": sequence.get("modified"),
            "sequence_version": audit.get("sequenceVersion"),
        },
        "entry_version": audit.get("entryVersion"),
        "last_annotation_update": audit.get("lastAnnotationUpdateDate"),
        "catalytic_status": "catalytic" if cfg["catalytic"] else "acatalytic",
        "catalytic_status_basis": {
            "function_texts": functions,
            "rhea_10748_assignment_present": reaction is not None,
        },
        "function_texts": functions,
        "subcellular_location_comments": [
            comment
            for comment in entry.get("comments", [])
            if comment.get("commentType") == "SUBCELLULAR LOCATION"
        ],
        "experimental_pdb_cross_reference_count": len(pdb_refs),
        "selected_structure_record_count": cfg["structures"],
        "chembl_target_ids": chembl_ids,
        "selected_activity_record_count": cfg["activities"],
        "selected_reaction_record_count": 1 if cfg["catalytic"] else 0,
        "canonical_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
        "provenance": [downloader.source_ref(rel_path)],
    }


def normalize_pdbe_entity(entity: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "molecule_type",
        "entity_id",
        "sample_preparation",
        "length",
        "number_of_copies",
        "in_chains",
        "in_struct_asyms",
        "mutation_flag",
        "weight",
        "ca_p_only",
        "synonym",
        "molecule_name",
        "gene_name",
        "source",
        "sequence",
        "pdb_sequence",
        "pdb_sequence_indices_with_multiple_residues",
    )
    return {key: entity.get(key) for key in keys}


def select_and_build_structures(
    configs: tuple[dict[str, Any], ...],
    uniprot_entries: dict[str, dict[str, Any]],
    downloader: SnapshotDownloader,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    mapping_payloads: dict[tuple[str, str], dict[str, Any]] = {}
    for cfg in configs:
        quota = cfg["structures"]
        if quota == 0:
            continue
        references = xrefs(uniprot_entries[cfg["accession"]], "PDB")
        if len(references) < quota:
            raise RuntimeError(f"Not enough UniProt PDB references for {cfg['gene']}: {len(references)} < {quota}")
        ordered = structure_candidate_order(references)
        gene_selected = 0
        for candidate in ordered:
            pdb_id = candidate["reference"]["id"].upper()
            rel_path = f"raw/pdbe/mappings/{pdb_id.lower()}.json"
            url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
            payload = downloader.json(rel_path, url)
            mapping = mapping_for_accession(payload, pdb_id, cfg["accession"])
            if not mapping or not mapping.get("mappings"):
                continue
            mapping_payloads[(cfg["accession"], pdb_id)] = payload
            selected.append(
                {
                    "cfg": cfg,
                    "candidate": candidate,
                    "pdb_id": pdb_id,
                    "rank_within_gene": gene_selected + 1,
                    "available_count": len(references),
                    "mapping_path": rel_path,
                }
            )
            gene_selected += 1
            if gene_selected == quota:
                break
        if gene_selected != quota:
            raise RuntimeError(f"Only {gene_selected}/{quota} exact SIFTS mappings were available for {cfg['gene']}")

    if len(selected) != EXPECTED_STRUCTURE_TOTAL:
        raise RuntimeError(f"Structure selection produced {len(selected)} records")

    pdb_ids = sorted({item["pdb_id"].lower() for item in selected})
    pdbe_indexes: dict[str, dict[str, tuple[Any, str]]] = {
        "summary": {},
        "experiment": {},
        "molecules": {},
    }
    chunk_size = 20
    for chunk_index, start in enumerate(range(0, len(pdb_ids), chunk_size), start=1):
        chunk = pdb_ids[start : start + chunk_size]
        joined = ",".join(chunk)
        for endpoint in ("summary", "experiment", "molecules"):
            rel_path = f"raw/pdbe/{endpoint}/chunk-{chunk_index:02d}.json"
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/{endpoint}/{joined}"
            payload = downloader.json(rel_path, url)
            for pdb_id, value in payload.items():
                pdbe_indexes[endpoint][pdb_id.lower()] = (value, rel_path)

    records: list[dict[str, Any]] = []
    for item in selected:
        cfg = item["cfg"]
        accession = cfg["accession"]
        pdb_id = item["pdb_id"]
        pdb_key = pdb_id.lower()
        mapping_payload = mapping_payloads[(accession, pdb_id)]
        mapped = mapping_for_accession(mapping_payload, pdb_id, accession)
        assert mapped is not None
        mappings = mapped["mappings"]
        summary_values, summary_path = pdbe_indexes["summary"][pdb_key]
        experiment_values, experiment_path = pdbe_indexes["experiment"][pdb_key]
        molecule_values, molecule_path = pdbe_indexes["molecules"][pdb_key]
        summary = summary_values[0] if summary_values else {}
        mapped_entities = {mapping.get("entity_id") for mapping in mappings}
        construct_entities = [
            normalize_pdbe_entity(entity)
            for entity in molecule_values
            if entity.get("entity_id") in mapped_entities
        ]
        if not construct_entities:
            raise RuntimeError(f"No PDBe molecule entity matched {pdb_id}/{accession}")
        props = property_map(item["candidate"]["reference"])
        provenance_paths = [
            f"raw/uniprot/entries/{accession}.json",
            item["mapping_path"],
            summary_path,
            experiment_path,
            molecule_path,
        ]
        records.append(
            {
                "record_id": f"STRUCTURE:PDB:{pdb_id}:UNIPROT:{accession}",
                "evidence_type": "experimental_structure_mapping",
                "family_id": FAMILY_ID,
                "protein": {
                    "gene_symbol": cfg["gene"],
                    "uniprot_accession": accession,
                    "taxon_id": 9606,
                },
                "structure": {
                    "pdb_id": pdb_id,
                    "title": summary.get("title"),
                    "experimental_method": summary.get("experimental_method"),
                    "experimental_method_class": summary.get("experimental_method_class"),
                    "deposition_date": summary.get("deposition_date"),
                    "release_date": summary.get("release_date"),
                    "revision_date": summary.get("revision_date"),
                    "entry_authors": summary.get("entry_authors"),
                    "number_of_entities": summary.get("number_of_entities"),
                    "assemblies": summary.get("assemblies"),
                    "experiments": experiment_values,
                },
                "uniprot_cross_reference": {
                    "method": props.get("Method"),
                    "resolution": props.get("Resolution"),
                    "chains": props.get("Chains"),
                },
                "sifts_mapping": {
                    "identifier": mapped.get("identifier"),
                    "name": mapped.get("name"),
                    "mappings": mappings,
                },
                "construct_entities": construct_entities,
                "selection": {
                    "strategy": "method-and-resolution-band-round-robin",
                    "method_bucket": item["candidate"]["method"],
                    "resolution_band": item["candidate"]["resolution_band"],
                    "reported_resolution_angstrom": item["candidate"]["resolution_angstrom"],
                    "rank_within_gene_slice": item["rank_within_gene"],
                    "gene_quota": cfg["structures"],
                    "available_uniprot_pdb_cross_references": item["available_count"],
                },
                "interpretation_boundary": {
                    "predicted_structure": False,
                    "co_crystal_ligand_implies_affinity": False,
                    "kinetics_or_flux_inferred": False,
                },
                "canonical_url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_key}",
                "provenance": [downloader.source_ref(path) for path in provenance_paths],
            }
        )
    return sorted(records, key=lambda record: record["record_id"])


def exact_target_component(target: dict[str, Any], accession: str) -> dict[str, Any]:
    if target.get("target_type") != "SINGLE PROTEIN":
        raise RuntimeError(f"ChEMBL target is not SINGLE PROTEIN: {target.get('target_chembl_id')}")
    if target.get("tax_id") != 9606 or target.get("organism") != "Homo sapiens":
        raise RuntimeError(f"ChEMBL target is not exact human target: {target.get('target_chembl_id')}")
    components = [
        component
        for component in target.get("target_components", [])
        if component.get("component_type") == "PROTEIN" and component.get("accession") == accession
    ]
    if len(components) != 1 or components[0].get("relationship") != "SINGLE PROTEIN":
        raise RuntimeError(f"ChEMBL component mapping mismatch for {accession}")
    protein_components = [
        component
        for component in target.get("target_components", [])
        if component.get("component_type") == "PROTEIN"
    ]
    if len(protein_components) != 1:
        raise RuntimeError(f"ChEMBL SINGLE PROTEIN target has {len(protein_components)} protein components")
    return components[0]


def activity_basics_valid(activity: dict[str, Any], target_id: str) -> bool:
    return all(
        (
            activity.get("target_chembl_id") == target_id,
            activity.get("activity_id") is not None,
            activity.get("assay_chembl_id"),
            activity.get("molecule_chembl_id"),
            activity.get("document_chembl_id"),
            activity.get("standard_type"),
            activity.get("standard_relation") is not None,
            activity.get("standard_value") is not None,
            activity.get("standard_units"),
            activity.get("potential_duplicate") in (0, False),
            activity.get("data_validity_comment") is None,
        )
    )


def assay_is_direct(assay: dict[str, Any], target_id: str) -> bool:
    return (
        assay.get("target_chembl_id") == target_id
        and assay.get("confidence_score") == 9
        and assay.get("relationship_type") == "D"
    )


def normalized_document(document: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "document_chembl_id",
        "chembl_release",
        "doc_type",
        "title",
        "authors",
        "journal",
        "journal_full_title",
        "year",
        "volume",
        "issue",
        "first_page",
        "last_page",
        "doi",
        "doi_chembl",
        "pubmed_id",
        "patent_id",
        "src_id",
    )
    return {key: document.get(key) for key in keys}


def normalized_molecule(molecule: dict[str, Any]) -> dict[str, Any]:
    structures = molecule.get("molecule_structures") or {}
    hierarchy = molecule.get("molecule_hierarchy") or {}
    return {
        "molecule_chembl_id": molecule.get("molecule_chembl_id"),
        "parent_molecule_chembl_id": hierarchy.get("parent_chembl_id"),
        "is_parent_form": hierarchy.get("parent_chembl_id") == molecule.get("molecule_chembl_id"),
        "pref_name": molecule.get("pref_name"),
        "molecule_type": molecule.get("molecule_type"),
        "structure_type": molecule.get("structure_type"),
        "chirality": molecule.get("chirality"),
        "standard_inchi": structures.get("standard_inchi"),
        "standard_inchi_key": structures.get("standard_inchi_key"),
        "canonical_smiles": structures.get("canonical_smiles"),
        "molecule_properties": molecule.get("molecule_properties"),
        "inorganic_flag": molecule.get("inorganic_flag"),
        "polymer_flag": molecule.get("polymer_flag"),
    }


def select_and_build_activities(
    configs: tuple[dict[str, Any], ...],
    uniprot_entries: dict[str, dict[str, Any]],
    downloader: SnapshotDownloader,
    chembl_release: str,
) -> list[dict[str, Any]]:
    active_configs = [cfg for cfg in configs if cfg["activities"]]
    target_specs: list[tuple[str, str]] = []
    targets_by_gene: dict[str, tuple[str, str]] = {}
    for cfg in active_configs:
        chembl_refs = xrefs(uniprot_entries[cfg["accession"]], "ChEMBL")
        if len(chembl_refs) != 1:
            raise RuntimeError(f"Expected one ChEMBL target for {cfg['gene']}, found {len(chembl_refs)}")
        target_id = chembl_refs[0]["id"]
        rel_path = f"raw/chembl/targets/{target_id}.json"
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json"
        targets_by_gene[cfg["gene"]] = (target_id, rel_path)
        target_specs.append((rel_path, url))
    fetched_targets = fetch_json_specs(downloader, target_specs)

    target_payloads: dict[str, dict[str, Any]] = {}
    target_components: dict[str, dict[str, Any]] = {}
    for cfg in active_configs:
        target_id, rel_path = targets_by_gene[cfg["gene"]]
        target = fetched_targets[rel_path]
        if target.get("target_chembl_id") != target_id:
            raise RuntimeError(f"ChEMBL target identifier mismatch: {target_id}")
        component = exact_target_component(target, cfg["accession"])
        target_payloads[cfg["gene"]] = target
        target_components[cfg["gene"]] = component

    # Start from ChEMBL's assay endpoint, where relationship/confidence are
    # first-class fields. This avoids treating a homolog-assigned activity as a
    # direct human-protein measurement merely because its target identifier was
    # normalized later. Within each exact target we inspect a deterministic
    # window of direct assays, then take the earliest eligible activity per
    # assay. No value or potency is used to rank candidates.
    direct_assays_by_gene: dict[str, list[dict[str, Any]]] = {}
    direct_assay_page_paths: dict[str, str] = {}
    for cfg in active_configs:
        target_id, _ = targets_by_gene[cfg["gene"]]
        query = urllib.parse.urlencode(
            {
                "target_chembl_id": target_id,
                "confidence_score": "9",
                "relationship_type": "D",
                "limit": "50",
                "order_by": "assay_chembl_id",
            }
        )
        rel_path = f"raw/chembl/direct-assays/{target_id}.json"
        url = f"https://www.ebi.ac.uk/chembl/api/data/assay.json?{query}"
        payload = downloader.json(rel_path, url)
        assays = [
            assay
            for assay in payload.get("assays", [])
            if assay_is_direct(assay, target_id)
        ]
        if len(assays) < cfg["activities"]:
            raise RuntimeError(f"Too few confidence-9 direct ChEMBL assays for {cfg['gene']}")
        direct_assays_by_gene[cfg["gene"]] = assays
        direct_assay_page_paths[cfg["gene"]] = rel_path

    assay_activity_specs: list[tuple[str, str]] = []
    inspected_assays_by_gene: dict[str, list[dict[str, Any]]] = {}
    for cfg in active_configs:
        # Ten fallbacks keep the normal source slice compact while covering
        # direct assays that contain no numeric standardized activity.
        inspected = direct_assays_by_gene[cfg["gene"]][: cfg["activities"] + 10]
        inspected_assays_by_gene[cfg["gene"]] = inspected
        for assay in inspected:
            assay_id = assay["assay_chembl_id"]
            target_id, _ = targets_by_gene[cfg["gene"]]
            query = urllib.parse.urlencode(
                {
                    "assay_chembl_id": assay_id,
                    "target_chembl_id": target_id,
                    "standard_value__isnull": "false",
                    "potential_duplicate": "0",
                    "data_validity_comment__isnull": "true",
                    "limit": "10",
                    "order_by": "activity_id",
                }
            )
            assay_activity_specs.append(
                (
                    f"raw/chembl/assay-activities/{assay_id}.json",
                    f"https://www.ebi.ac.uk/chembl/api/data/activity.json?{query}",
                )
            )
    assay_activity_pages = fetch_json_specs(downloader, assay_activity_specs)

    detail_candidates: dict[str, list[dict[str, Any]]] = {}
    detail_specs: list[tuple[str, str]] = []
    for cfg in active_configs:
        target_id, _ = targets_by_gene[cfg["gene"]]
        valid: list[dict[str, Any]] = []
        for assay in inspected_assays_by_gene[cfg["gene"]]:
            assay_id = assay["assay_chembl_id"]
            candidate_path = f"raw/chembl/assay-activities/{assay_id}.json"
            page = assay_activity_pages[candidate_path]
            eligible = [
                activity
                for activity in page.get("activities", [])
                if activity_basics_valid(activity, target_id)
                and activity.get("assay_chembl_id") == assay_id
            ]
            if eligible:
                valid.append(
                    {
                        "activity": sorted(eligible, key=lambda item: int(item["activity_id"]))[0],
                        "assay": assay,
                        "candidate_path": candidate_path,
                    }
                )
            if len(valid) >= cfg["activities"] + 2:
                break
        if len(valid) < cfg["activities"]:
            raise RuntimeError(f"Too few strictly admissible direct-assay activities for {cfg['gene']}")
        detail_candidates[cfg["gene"]] = valid
        for candidate in valid:
            activity = candidate["activity"]
            activity_id = activity["activity_id"]
            assay_id = activity["assay_chembl_id"]
            molecule_id = activity["molecule_chembl_id"]
            document_id = activity["document_chembl_id"]
            detail_specs.extend(
                (
                    (
                        f"raw/chembl/activities/{activity_id}.json",
                        f"https://www.ebi.ac.uk/chembl/api/data/activity/{activity_id}.json",
                    ),
                    (
                        f"raw/chembl/assays/{assay_id}.json",
                        f"https://www.ebi.ac.uk/chembl/api/data/assay/{assay_id}.json",
                    ),
                    (
                        f"raw/chembl/molecules/{molecule_id}.json",
                        f"https://www.ebi.ac.uk/chembl/api/data/molecule/{molecule_id}.json",
                    ),
                    (
                        f"raw/chembl/documents/{document_id}.json",
                        f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json",
                    ),
                )
            )
    details = fetch_json_specs(downloader, detail_specs)

    records: list[dict[str, Any]] = []
    for cfg in active_configs:
        target_id, target_path = targets_by_gene[cfg["gene"]]
        target = target_payloads[cfg["gene"]]
        component = target_components[cfg["gene"]]
        chosen: list[
            tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]
        ] = []
        for candidate in detail_candidates[cfg["gene"]]:
            candidate_activity = candidate["activity"]
            activity_id = candidate_activity["activity_id"]
            molecule_id = candidate_activity["molecule_chembl_id"]
            document_id = candidate_activity["document_chembl_id"]
            activity = details[f"raw/chembl/activities/{activity_id}.json"]
            assay = details[f"raw/chembl/assays/{candidate_activity['assay_chembl_id']}.json"]
            molecule = details[f"raw/chembl/molecules/{molecule_id}.json"]
            document = details[f"raw/chembl/documents/{document_id}.json"]
            if not activity_basics_valid(activity, target_id):
                continue
            if not assay_is_direct(assay, target_id):
                continue
            structures = molecule.get("molecule_structures") or {}
            if molecule.get("molecule_chembl_id") != molecule_id or not structures.get("standard_inchi_key"):
                continue
            if document.get("document_chembl_id") != document_id:
                continue
            chosen.append((activity, assay, molecule, document, candidate))
            if len(chosen) == cfg["activities"]:
                break
        if len(chosen) != cfg["activities"]:
            raise RuntimeError(f"Could not select exact compound-identified activity quota for {cfg['gene']}")

        for rank, (activity, assay, molecule, document, candidate) in enumerate(chosen, start=1):
            activity_id = activity["activity_id"]
            assay_id = activity["assay_chembl_id"]
            molecule_id = activity["molecule_chembl_id"]
            document_id = activity["document_chembl_id"]
            assay_path = f"raw/chembl/assays/{assay_id}.json"
            source_paths = [
                f"raw/uniprot/entries/{cfg['accession']}.json",
                target_path,
                direct_assay_page_paths[cfg["gene"]],
                candidate["candidate_path"],
                f"raw/chembl/activities/{activity_id}.json",
                assay_path,
                f"raw/chembl/molecules/{molecule_id}.json",
                f"raw/chembl/documents/{document_id}.json",
            ]
            comparison_refusals = []
            if not cfg["catalytic"]:
                comparison_refusals.append("CATALYTIC_STATUS_MISMATCH")
            if not assay.get("assay_parameters"):
                comparison_refusals.append("CONDITION_MISSING")
            if activity.get("standard_relation") != "=":
                comparison_refusals.append("CENSORED_RELATION")
            if activity.get("standard_units") not in {"M", "mM", "uM", "nM", "pM", "fM"}:
                comparison_refusals.append("UNSUPPORTED_UNIT")
            records.append(
                {
                    "record_id": f"ACTIVITY:{chembl_release}:{activity_id}",
                    "evidence_type": "assay_activity",
                    "family_id": FAMILY_ID,
                    "protein": {
                        "gene_symbol": cfg["gene"],
                        "uniprot_accession": cfg["accession"],
                        "catalytic_status": "catalytic" if cfg["catalytic"] else "acatalytic",
                    },
                    "target": {
                        "target_chembl_id": target_id,
                        "target_type": target.get("target_type"),
                        "pref_name": target.get("pref_name"),
                        "organism": target.get("organism"),
                        "tax_id": target.get("tax_id"),
                        "component": {
                            "component_id": component.get("component_id"),
                            "component_type": component.get("component_type"),
                            "accession": component.get("accession"),
                            "relationship": component.get("relationship"),
                            "component_description": component.get("component_description"),
                        },
                    },
                    "assay": assay,
                    "measurement": {
                        "reported": {
                            "type": activity.get("type"),
                            "relation": activity.get("relation"),
                            "value": activity.get("value"),
                            "upper_value": activity.get("upper_value"),
                            "units": activity.get("units"),
                            "text_value": activity.get("text_value"),
                        },
                        "standardized_by_source": {
                            "type": activity.get("standard_type"),
                            "relation": activity.get("standard_relation"),
                            "value": activity.get("standard_value"),
                            "upper_value": activity.get("standard_upper_value"),
                            "units": activity.get("standard_units"),
                            "text_value": activity.get("standard_text_value"),
                            "flag": activity.get("standard_flag"),
                        },
                        "pchembl_value_from_source": activity.get("pchembl_value"),
                        "bao_endpoint": activity.get("bao_endpoint"),
                        "bao_format": activity.get("bao_format"),
                        "bao_label": activity.get("bao_label"),
                        "activity_properties": activity.get("activity_properties"),
                        "assay_variant_accession": activity.get("assay_variant_accession"),
                        "assay_variant_mutation": activity.get("assay_variant_mutation"),
                        "data_validity_comment": activity.get("data_validity_comment"),
                        "data_validity_description": activity.get("data_validity_description"),
                        "potential_duplicate": activity.get("potential_duplicate"),
                    },
                    "compound": normalized_molecule(molecule),
                    "document": normalized_document(document),
                    "selection": {
                        "strategy": "ascending-direct-assay-id-then-activity-id",
                        "rank_within_target_slice": rank,
                        "target_quota": cfg["activities"],
                        "candidate_activity_id": candidate["activity"]["activity_id"],
                        "excluded_source_flags": ["potential_duplicate", "data_validity_comment"],
                    },
                    "comparability": {
                        "prequalified": False,
                        "record_local_refusal_codes": comparison_refusals,
                        "note": "No cross-assay comparison is authorized by this record. Evaluate target, construct, endpoint, relation, units, assay format, conditions, and exact molecular form first.",
                    },
                    "canonical_urls": {
                        "activity_api": f"https://www.ebi.ac.uk/chembl/api/data/activity/{activity_id}.json",
                        "assay_api": f"https://www.ebi.ac.uk/chembl/api/data/assay/{assay_id}.json",
                        "target_api": f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json",
                        "molecule_api": f"https://www.ebi.ac.uk/chembl/api/data/molecule/{molecule_id}.json",
                        "document_api": f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json",
                    },
                    "provenance": [downloader.source_ref(path) for path in source_paths],
                }
            )
    if len(records) != EXPECTED_ACTIVITY_TOTAL:
        raise RuntimeError(f"Activity selection produced {len(records)} records")
    return sorted(records, key=lambda record: record["record_id"])


def build_reactions(
    configs: tuple[dict[str, Any], ...],
    uniprot_entries: dict[str, dict[str, Any]],
    rhea_row: dict[str, str],
    downloader: SnapshotDownloader,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for cfg in configs:
        if not cfg["catalytic"]:
            continue
        reaction = reaction_with_rhea(uniprot_entries[cfg["accession"]], RHEA_MASTER_ID)
        if reaction is None:
            raise RuntimeError(f"Missing UniProt/Rhea assignment for {cfg['gene']}")
        source_paths = [
            f"raw/uniprot/entries/{cfg['accession']}.json",
            "raw/rhea/reaction-10748.tsv",
            "raw/rhea/rhea-release.properties",
        ]
        records.append(
            {
                "record_id": f"REACTION:UNIPROT:{cfg['accession']}:RHEA:10748",
                "evidence_type": "reaction_assignment",
                "family_id": FAMILY_ID,
                "protein": {
                    "gene_symbol": cfg["gene"],
                    "uniprot_accession": cfg["accession"],
                    "taxon_id": 9606,
                    "catalytic_status": "catalytic",
                },
                "uniprot_assignment": {
                    "reaction_name": reaction.get("name"),
                    "ec_number": reaction.get("ecNumber"),
                    "cross_references": reaction.get("reactionCrossReferences"),
                    "evidences": reaction.get("evidences"),
                },
                "rhea_reaction": {
                    "rhea_id": rhea_row.get("Reaction identifier"),
                    "equation": rhea_row.get("Equation"),
                    "enzyme_class": rhea_row.get("Enzyme class"),
                    "linked_enzyme_count_from_snapshot": int(rhea_row["Enzymes"])
                    if rhea_row.get("Enzymes", "").isdigit()
                    else None,
                },
                "interpretation_boundary": {
                    "directional_flux_claim": False,
                    "kinetic_parameter_claim": False,
                    "tissue_activity_claim": False,
                },
                "canonical_urls": {
                    "uniprot": f"https://www.uniprot.org/uniprotkb/{cfg['accession']}/entry",
                    "rhea": "https://www.rhea-db.org/rhea/10748",
                },
                "provenance": [downloader.source_ref(path) for path in source_paths],
            }
        )
    if len(records) != EXPECTED_REACTION_TOTAL:
        raise RuntimeError(f"Reaction build produced {len(records)} records")
    return sorted(records, key=lambda record: record["record_id"])


def source_licenses(downloader: SnapshotDownloader, chembl_release: str) -> dict[str, Any]:
    return {
        "aggregate_license": None,
        "aggregate_license_note": "No blanket license is asserted for this mixed-source aggregation. Reusers must comply with each source license and preserve record provenance.",
        "generator_code_license": "Apache-2.0",
        "record_licenses": {
            record_type: {"license": license_id, "license_url": LICENSE_URLS[license_id]}
            for record_type, license_id in RECORD_LICENSES.items()
        },
        "sources": [
            {
                "source": "UniProtKB",
                "license": "CC-BY-4.0",
                "official_license_url": LICENSE_URLS["CC-BY-4.0"],
                "license_text": "LICENSES/CC-BY-4.0.txt",
                "terms_url": "https://www.uniprot.org/help/license",
                "snapshot": downloader.source_ref("raw/licenses/uniprot-license.txt"),
                "attribution_required": True,
            },
            {
                "source": "Rhea",
                "license": "CC-BY-4.0",
                "official_license_url": LICENSE_URLS["CC-BY-4.0"],
                "license_text": "LICENSES/CC-BY-4.0.txt",
                "terms_url": "https://ftp.expasy.org/databases/rhea/LICENSE.txt",
                "snapshot": downloader.source_ref("raw/licenses/rhea-license.txt"),
                "attribution_required": True,
            },
            {
                "source": "ChEMBL",
                "license": "CC-BY-SA-3.0",
                "official_license_url": LICENSE_URLS["CC-BY-SA-3.0"],
                "license_text": "LICENSES/CC-BY-SA-3.0.txt",
                "terms_url": "https://chembl.gitbook.io/chembl-interface-documentation/about",
                "snapshot": downloader.source_ref("raw/licenses/chembl-license-extract.txt"),
                "attribution_required": True,
                "share_alike": True,
                "citation": chembl_citation(chembl_release),
                "applies_to": [
                    "corpus/activities.json",
                    "assay_activity records in corpus/corpus.json",
                    "raw/chembl/**",
                ],
            },
            {
                "source": "PDBe and wwPDB public data",
                "license": "CC0-1.0",
                "official_license_url": LICENSE_URLS["CC0-1.0"],
                "license_text": "LICENSES/CC0-1.0.txt",
                "terms_url": "https://www.ebi.ac.uk/pdbe/about/public-data-access-statement",
                "snapshot": downloader.source_ref("raw/licenses/pdbe-license-extract.txt"),
                "attribution_requested_as_scientific_practice": True,
            },
        ],
    }


def with_record_license(record: dict[str, Any]) -> dict[str, Any]:
    license_id = RECORD_LICENSES[record["evidence_type"]]
    return {**record, "license": license_id, "license_url": LICENSE_URLS[license_id]}


def write_checksums(directory: Path, output_name: str = "checksums.sha256") -> None:
    rows = []
    for path in sorted(item for item in directory.rglob("*") if item.is_file() and item.name != output_name):
        rows.append(f"{sha256_bytes(path.read_bytes())}  {path.relative_to(directory).as_posix()}")
    write_text(directory / output_name, "\n".join(rows) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--refresh", action="store_true", help="replace every cached HTTP snapshot")
    mode.add_argument("--offline", action="store_true", help="forbid network and require a complete cache")
    args = parser.parse_args()

    if sum(cfg["structures"] for cfg in CORE_PROTEINS) != EXPECTED_STRUCTURE_TOTAL:
        raise RuntimeError("Pinned structure quotas do not sum to 80")
    if sum(cfg["activities"] for cfg in CORE_PROTEINS) != EXPECTED_ACTIVITY_TOTAL:
        raise RuntimeError("Pinned activity quotas do not sum to 28")
    if sum(1 for cfg in CORE_PROTEINS if cfg["catalytic"]) != EXPECTED_REACTION_TOTAL:
        raise RuntimeError("Pinned catalytic set does not sum to 12")

    downloader = SnapshotDownloader(refresh=args.refresh, offline=args.offline)

    for rel_path, url, accept in LICENSE_SOURCES:
        downloader.fetch(rel_path, url, accept=accept)
    for rel_path, url, phrase in LICENSE_PAGE_EXTRACTS:
        downloader.extract(rel_path, url, phrase)

    chembl_status = downloader.json(
        "raw/chembl/status.json", "https://www.ebi.ac.uk/chembl/api/data/status.json"
    )
    chembl_release = chembl_status.get("chembl_db_version")
    if not chembl_release:
        raise RuntimeError("ChEMBL status did not expose a release identifier")

    rhea_properties_text = downloader.text(
        "raw/rhea/rhea-release.properties",
        "https://ftp.expasy.org/databases/rhea/rhea-release.properties",
    )
    rhea_properties = {}
    for line in rhea_properties_text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            rhea_properties[key.strip()] = value.strip()
    rhea_query = urllib.parse.urlencode({"query": "rhea:10748", "format": "tsv"})
    rhea_text = downloader.text(
        "raw/rhea/reaction-10748.tsv",
        f"https://www.rhea-db.org/rhea?{rhea_query}",
        accept="text/tab-separated-values",
    )
    rhea_row = parse_rhea_tsv(rhea_text)

    search_terms = " OR ".join(f"gene_exact:{cfg['gene']}" for cfg in CORE_PROTEINS)
    search_query = f"organism_id:9606 AND reviewed:true AND ({search_terms})"
    family_search_url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": search_query, "format": "json", "size": "50"}
    )
    family_search = downloader.json("raw/uniprot/family-search.json", family_search_url)
    primary_gene_accessions = {
        primary_gene(entry): entry.get("primaryAccession") for entry in family_search.get("results", [])
    }
    for cfg in CORE_PROTEINS:
        if primary_gene_accessions.get(cfg["gene"]) != cfg["accession"]:
            raise RuntimeError(f"UniProt search did not resolve pinned primary identity for {cfg['gene']}")

    uniprot_specs = [
        (
            f"raw/uniprot/entries/{cfg['accession']}.json",
            f"https://rest.uniprot.org/uniprotkb/{cfg['accession']}.json",
        )
        for cfg in CORE_PROTEINS
    ]
    fetched_uniprot = fetch_json_specs(downloader, uniprot_specs)
    uniprot_entries = {
        cfg["accession"]: fetched_uniprot[f"raw/uniprot/entries/{cfg['accession']}.json"]
        for cfg in CORE_PROTEINS
    }
    proteins = [normalize_protein(cfg, uniprot_entries[cfg["accession"]], downloader) for cfg in CORE_PROTEINS]
    downloader.save_manifest()

    structures = select_and_build_structures(CORE_PROTEINS, uniprot_entries, downloader)
    downloader.save_manifest()
    activities = select_and_build_activities(
        CORE_PROTEINS, uniprot_entries, downloader, str(chembl_release)
    )
    downloader.save_manifest()
    reactions = build_reactions(CORE_PROTEINS, uniprot_entries, rhea_row, downloader)
    structures = [with_record_license(record) for record in structures]
    activities = [with_record_license(record) for record in activities]
    reactions = [with_record_license(record) for record in reactions]

    downloader.save_manifest()
    write_checksums(ROOT / "raw")

    uniprot_release_values = {
        entry.get("response_headers", {}).get("x-uniprot-release")
        for path, entry in downloader.entries.items()
        if path.startswith("raw/uniprot/entries/")
    }
    uniprot_release_values.discard(None)
    if len(uniprot_release_values) != 1:
        raise RuntimeError(f"UniProt snapshots span releases: {sorted(uniprot_release_values)}")
    uniprot_release = next(iter(uniprot_release_values))
    retrieval_times = sorted(entry["retrieved_at"] for entry in downloader.entries.values())
    generated_at = retrieval_times[-1]
    licenses = source_licenses(downloader, str(chembl_release))
    all_records = structures + activities + reactions
    release_digest = sha256_bytes(canonical_bytes({"proteins": proteins, "evidence_records": all_records}))[:16]
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "corpus_id": FAMILY_ID + "-evidence-v1",
        "release_id": f"haca-evidence-v1-{release_digest}",
        "title": "Human alpha-carbonic anhydrase evidence corpus",
        "generated_from_snapshots_at": generated_at,
        "record_counts": {
            "proteins": len(proteins),
            "experimental_structure_mapping": len(structures),
            "assay_activity": len(activities),
            "reaction_assignment": len(reactions),
            "total_evidence_records": len(all_records),
        },
        "source_releases": {
            "uniprot": uniprot_release,
            "chembl": chembl_release,
            "chembl_release_date": chembl_status.get("chembl_release_date"),
            "rhea": rhea_properties.get("rhea.release.number"),
            "rhea_release_date": rhea_properties.get("rhea.release.date"),
            "pdbe": "live-api-snapshot; per-entry revision_date retained",
        },
        "scope": {
            "organism": "Homo sapiens",
            "taxon_id": 9606,
            "core_paralog_count": 15,
            "catalytic_paralog_count": 12,
            "acatalytic_related_paralog_count": 3,
            "clinical_use": False,
            "kinetic_or_flux_model": False,
        },
        "selection_policy_document": "README.md#selection-policy",
        "data_dictionary": "docs/data-dictionary.md",
    }

    family_file = {"metadata": metadata, "proteins": proteins}
    structures_file = {
        "schema_version": SCHEMA_VERSION,
        "evidence_type": "experimental_structure_mapping",
        "count": len(structures),
        "records": structures,
    }
    activities_file = {
        "schema_version": SCHEMA_VERSION,
        "evidence_type": "assay_activity",
        "license": RECORD_LICENSES["assay_activity"],
        "license_url": LICENSE_URLS[RECORD_LICENSES["assay_activity"]],
        "citation": chembl_citation(str(chembl_release)),
        "count": len(activities),
        "records": activities,
    }
    reactions_file = {
        "schema_version": SCHEMA_VERSION,
        "evidence_type": "reaction_assignment",
        "count": len(reactions),
        "records": reactions,
    }
    corpus_file = {
        "metadata": metadata,
        "proteins": proteins,
        "evidence_records": all_records,
        "licenses": licenses,
    }
    source_snapshots = {
        "schema_version": SCHEMA_VERSION,
        "source_releases": metadata["source_releases"],
        "retrieval_manifest": "raw/retrieval-manifest.json",
        "raw_checksums": "raw/checksums.sha256",
        "snapshots": [downloader.entries[path] for path in sorted(downloader.entries)],
    }

    write_json(ROOT / "corpus/family.json", family_file)
    write_json(ROOT / "corpus/structures.json", structures_file)
    write_json(ROOT / "corpus/activities.json", activities_file)
    write_json(ROOT / "corpus/reactions.json", reactions_file)
    write_json(ROOT / "corpus/licenses.json", licenses)
    write_json(ROOT / "corpus/source-snapshots.json", source_snapshots)
    write_json(ROOT / "corpus/corpus.json", corpus_file)
    write_checksums(ROOT / "corpus")

    print(
        json.dumps(
            {
                "release_id": metadata["release_id"],
                "proteins": len(proteins),
                "structures": len(structures),
                "activities": len(activities),
                "reactions": len(reactions),
                "evidence_total": len(all_records),
                "uniprot_release": uniprot_release,
                "chembl_release": chembl_release,
                "rhea_release": rhea_properties.get("rhea.release.number"),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
