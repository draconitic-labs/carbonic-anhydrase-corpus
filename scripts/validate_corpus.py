#!/usr/bin/env python3
"""Validate scientific identity, counts, provenance, and release integrity."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "CA1": {"accession": "P00915", "catalytic": True, "structures": 10, "activities": 2},
    "CA2": {"accession": "P00918", "catalytic": True, "structures": 20, "activities": 4},
    "CA3": {"accession": "P07451", "catalytic": True, "structures": 5, "activities": 2},
    "CA4": {"accession": "P22748", "catalytic": True, "structures": 10, "activities": 2},
    "CA5A": {"accession": "P35218", "catalytic": True, "structures": 0, "activities": 2},
    "CA5B": {"accession": "Q9Y2D0", "catalytic": True, "structures": 0, "activities": 2},
    "CA6": {"accession": "P23280", "catalytic": True, "structures": 1, "activities": 2},
    "CA7": {"accession": "P43166", "catalytic": True, "structures": 8, "activities": 2},
    "CA8": {"accession": "P35219", "catalytic": False, "structures": 1, "activities": 0},
    "CA9": {"accession": "Q16790", "catalytic": True, "structures": 8, "activities": 2},
    "CA10": {"accession": "Q9NS85", "catalytic": False, "structures": 0, "activities": 0},
    "CA11": {"accession": "O75493", "catalytic": False, "structures": 0, "activities": 2},
    "CA12": {"accession": "O43570", "catalytic": True, "structures": 8, "activities": 2},
    "CA13": {"accession": "Q8N1Q1", "catalytic": True, "structures": 7, "activities": 2},
    "CA14": {"accession": "Q9ULX7", "catalytic": True, "structures": 2, "activities": 2},
}
SUPPORTED_COMPARISON_UNITS = {"M", "mM", "uM", "nM", "pM", "fM"}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def validate_checksum_file(directory: Path) -> None:
    checksum_path = directory / "checksums.sha256"
    require(checksum_path.is_file(), f"Missing {checksum_path}")
    seen: set[str] = set()
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        expected_hash, rel_path = line.split("  ", 1)
        path = directory / rel_path
        require(path.is_file(), f"Checksummed file is missing: {path}")
        require(sha256_file(path) == expected_hash, f"Checksum mismatch: {path}")
        seen.add(rel_path)
    actual = {
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file() and path.name != "checksums.sha256"
    }
    require(seen == actual, f"Checksum coverage mismatch under {directory}: {sorted(actual ^ seen)}")


def validate_provenance(records: list[dict[str, Any]], snapshot_index: dict[str, dict[str, Any]]) -> None:
    for record in records:
        provenance = record.get("provenance") or []
        require(provenance, f"Record lacks provenance: {record.get('record_id')}")
        for source in provenance:
            path_text = source.get("snapshot_path")
            require(path_text in snapshot_index, f"Unknown source snapshot {path_text}")
            manifest = snapshot_index[path_text]
            require(source.get("snapshot_sha256") == manifest["sha256"], f"Source hash drift: {path_text}")
            require(source.get("url") == manifest["url"], f"Source URL drift: {path_text}")
            require(source.get("retrieved_at") == manifest["retrieved_at"], f"Retrieval time drift: {path_text}")
            path = ROOT / path_text
            require(path.is_file(), f"Source snapshot missing: {path_text}")
            require(sha256_file(path) == manifest["sha256"], f"Source snapshot changed: {path_text}")


def validate_proteins(proteins: list[dict[str, Any]]) -> None:
    require(len(proteins) == 15, "Expected 15 protein dossiers")
    require(len({item["protein_id"] for item in proteins}) == 15, "Protein IDs are not unique")
    by_gene = {item["gene_symbol"]: item for item in proteins}
    require(set(by_gene) == set(EXPECTED), "Core gene set mismatch")
    for gene, expected in EXPECTED.items():
        protein = by_gene[gene]
        require(protein["uniprot_accession"] == expected["accession"], f"Accession mismatch: {gene}")
        require(protein["taxon_id"] == 9606 and protein["organism"] == "Homo sapiens", f"Taxon mismatch: {gene}")
        require(protein["review_status"] == "reviewed", f"Unreviewed entry: {gene}")
        require(protein["sequence"]["value"], f"Missing sequence: {gene}")
        require(protein["sequence"]["crc64"], f"Missing sequence CRC64: {gene}")
        require(protein["sequence"]["sequence_version"] is not None, f"Missing sequence version: {gene}")
        expected_status = "catalytic" if expected["catalytic"] else "acatalytic"
        require(protein["catalytic_status"] == expected_status, f"Catalytic status mismatch: {gene}")
        require(
            protein["catalytic_status_basis"]["rhea_10748_assignment_present"] == expected["catalytic"],
            f"Rhea status mismatch: {gene}",
        )
        require(protein["selected_structure_record_count"] == expected["structures"], f"Structure quota mismatch: {gene}")
        require(protein["selected_activity_record_count"] == expected["activities"], f"Activity quota mismatch: {gene}")
        if expected["structures"] == 0:
            require(protein["experimental_pdb_cross_reference_count"] == 0, f"Hidden structures for {gene}")


def validate_structures(records: list[dict[str, Any]]) -> None:
    require(len(records) == 80, "Expected 80 structure records")
    counts = Counter(record["protein"]["gene_symbol"] for record in records)
    for gene, expected in EXPECTED.items():
        require(counts[gene] == expected["structures"], f"Structure count mismatch: {gene}")
    for record in records:
        require(record["evidence_type"] == "experimental_structure_mapping", "Structure type mismatch")
        gene = record["protein"]["gene_symbol"]
        accession = EXPECTED[gene]["accession"]
        require(record["protein"]["uniprot_accession"] == accession, f"Structure accession mismatch: {gene}")
        require(record["protein"]["taxon_id"] == 9606, f"Structure taxon mismatch: {gene}")
        require(record["structure"]["pdb_id"], f"Missing PDB ID: {record['record_id']}")
        require(record["structure"]["experimental_method"], f"Missing experimental method: {record['record_id']}")
        require(record["interpretation_boundary"]["predicted_structure"] is False, "Predicted structure admitted")
        mappings = record["sifts_mapping"]["mappings"]
        require(mappings, f"No SIFTS mapping: {record['record_id']}")
        entity_ids = {mapping["entity_id"] for mapping in mappings}
        require(all(mapping.get("chain_id") for mapping in mappings), f"Missing mapped chain: {record['record_id']}")
        require(all(mapping.get("struct_asym_id") for mapping in mappings), f"Missing asym ID: {record['record_id']}")
        require(all(mapping.get("identity") is not None for mapping in mappings), f"Missing identity: {record['record_id']}")
        require(all(mapping.get("coverage") is not None for mapping in mappings), f"Missing coverage: {record['record_id']}")
        construct_ids = {entity["entity_id"] for entity in record["construct_entities"]}
        require(entity_ids <= construct_ids, f"Mapped construct entity is absent: {record['record_id']}")
        require(record["selection"]["strategy"] == "method-and-resolution-band-round-robin", "Selection strategy drift")


def validate_activities(records: list[dict[str, Any]]) -> None:
    require(len(records) == 28, "Expected 28 activity records")
    counts = Counter(record["protein"]["gene_symbol"] for record in records)
    for gene, expected in EXPECTED.items():
        require(counts[gene] == expected["activities"], f"Activity count mismatch: {gene}")
    for record in records:
        require(record["evidence_type"] == "assay_activity", "Activity type mismatch")
        gene = record["protein"]["gene_symbol"]
        accession = EXPECTED[gene]["accession"]
        target = record["target"]
        assay = record["assay"]
        measurement = record["measurement"]
        compound = record["compound"]
        require(target["target_type"] == "SINGLE PROTEIN", f"Target scope mismatch: {record['record_id']}")
        require(target["tax_id"] == 9606 and target["organism"] == "Homo sapiens", f"Target taxon mismatch: {record['record_id']}")
        require(target["component"]["accession"] == accession, f"Target component mismatch: {record['record_id']}")
        require(target["component"]["relationship"] == "SINGLE PROTEIN", f"Component relationship mismatch: {record['record_id']}")
        require(assay["target_chembl_id"] == target["target_chembl_id"], f"Assay target mismatch: {record['record_id']}")
        require(assay["confidence_score"] == 9, f"Non-direct confidence admitted: {record['record_id']}")
        require(assay["relationship_type"] == "D", f"Non-direct relationship admitted: {record['record_id']}")
        standard = measurement["standardized_by_source"]
        require(standard["type"], f"Missing standard type: {record['record_id']}")
        require(standard["relation"] is not None, f"Missing relation: {record['record_id']}")
        require(standard["value"] is not None, f"Missing value: {record['record_id']}")
        require(standard["units"], f"Missing units: {record['record_id']}")
        require(measurement["potential_duplicate"] in (0, False), f"Potential duplicate admitted: {record['record_id']}")
        require(measurement["data_validity_comment"] is None, f"Validity-flagged activity admitted: {record['record_id']}")
        require(compound["molecule_chembl_id"], f"Missing molecule ID: {record['record_id']}")
        require(compound["parent_molecule_chembl_id"], f"Missing parent molecule ID: {record['record_id']}")
        require(compound["standard_inchi_key"], f"Missing InChIKey: {record['record_id']}")
        require(record["comparability"]["prequalified"] is False, "Activity incorrectly prequalified")
        refusal_codes = set(record["comparability"]["record_local_refusal_codes"])
        require("CONDITION_MISSING" in refusal_codes, f"Missing condition refusal: {record['record_id']}")
        if not EXPECTED[gene]["catalytic"]:
            require("CATALYTIC_STATUS_MISMATCH" in refusal_codes, f"Missing acatalytic refusal: {record['record_id']}")
        if standard["relation"] != "=":
            require("CENSORED_RELATION" in refusal_codes, f"Missing censored-relation refusal: {record['record_id']}")
        if standard["units"] not in SUPPORTED_COMPARISON_UNITS:
            require("UNSUPPORTED_UNIT" in refusal_codes, f"Missing unit refusal: {record['record_id']}")

        raw_path = next(
            source["snapshot_path"]
            for source in record["provenance"]
            if source["snapshot_path"].startswith("raw/chembl/activities/")
        )
        raw = load_json(ROOT / raw_path)
        require(str(raw["activity_id"]) in record["record_id"], f"Raw activity ID mismatch: {record['record_id']}")
        require(raw.get("standard_relation") == standard["relation"], f"Relation not preserved: {record['record_id']}")
        require(raw.get("standard_value") == standard["value"], f"Value not preserved: {record['record_id']}")
        require(raw.get("standard_units") == standard["units"], f"Units not preserved: {record['record_id']}")
        require(raw.get("assay_variant_accession") == measurement["assay_variant_accession"], f"Variant accession drift: {record['record_id']}")
        require(raw.get("assay_variant_mutation") == measurement["assay_variant_mutation"], f"Variant mutation drift: {record['record_id']}")


def validate_reactions(records: list[dict[str, Any]]) -> None:
    require(len(records) == 12, "Expected 12 reaction assignments")
    genes = {record["protein"]["gene_symbol"] for record in records}
    require(genes == {gene for gene, value in EXPECTED.items() if value["catalytic"]}, "Reaction gene set mismatch")
    for record in records:
        gene = record["protein"]["gene_symbol"]
        require(record["evidence_type"] == "reaction_assignment", "Reaction type mismatch")
        require(record["protein"]["uniprot_accession"] == EXPECTED[gene]["accession"], f"Reaction accession mismatch: {gene}")
        require(record["rhea_reaction"]["rhea_id"] == "RHEA:10748", f"Rhea ID mismatch: {gene}")
        require(record["rhea_reaction"]["equation"] == "hydrogencarbonate + H(+) = CO2 + H2O", f"Equation mismatch: {gene}")
        refs = record["uniprot_assignment"]["cross_references"]
        require(any(ref.get("database") == "Rhea" and ref.get("id") == "RHEA:10748" for ref in refs), f"UniProt Rhea evidence missing: {gene}")


def main() -> None:
    corpus = load_json(ROOT / "corpus/corpus.json")
    family = load_json(ROOT / "corpus/family.json")
    structures_file = load_json(ROOT / "corpus/structures.json")
    activities_file = load_json(ROOT / "corpus/activities.json")
    reactions_file = load_json(ROOT / "corpus/reactions.json")
    snapshots = load_json(ROOT / "corpus/source-snapshots.json")
    raw_manifest = load_json(ROOT / "raw/retrieval-manifest.json")

    metadata = corpus["metadata"]
    require(metadata["schema_version"] == "1.0.0", "Unexpected schema version")
    require(metadata["corpus_id"] == "human-alpha-carbonic-anhydrase-evidence-v1", "Corpus ID drift")
    require(metadata["record_counts"] == {
        "proteins": 15,
        "experimental_structure_mapping": 80,
        "assay_activity": 28,
        "reaction_assignment": 12,
        "total_evidence_records": 120,
    }, "Metadata counts mismatch")
    require(metadata["scope"]["clinical_use"] is False, "Clinical-use boundary changed")
    require(metadata["scope"]["kinetic_or_flux_model"] is False, "Kinetics boundary changed")
    require(metadata["source_releases"]["uniprot"] == "2026_02", "Unexpected UniProt release")
    require(metadata["source_releases"]["chembl"] == "ChEMBL_37", "Unexpected ChEMBL release")
    require(metadata["source_releases"]["rhea"] == "141", "Unexpected Rhea release")

    proteins = corpus["proteins"]
    records = corpus["evidence_records"]
    require(len(records) == 120, "Corpus evidence total mismatch")
    require(len({record["record_id"] for record in records}) == 120, "Evidence IDs are not unique")
    structures = [record for record in records if record["evidence_type"] == "experimental_structure_mapping"]
    activities = [record for record in records if record["evidence_type"] == "assay_activity"]
    reactions = [record for record in records if record["evidence_type"] == "reaction_assignment"]
    require(family["proteins"] == proteins, "family.json is not identical to corpus proteins")
    require(structures_file["records"] == structures, "structures.json diverges from corpus")
    require(activities_file["records"] == activities, "activities.json diverges from corpus")
    require(reactions_file["records"] == reactions, "reactions.json diverges from corpus")

    validate_proteins(proteins)
    validate_structures(structures)
    validate_activities(activities)
    validate_reactions(reactions)

    snapshot_entries = raw_manifest["entries"]
    snapshot_index = {entry["path"]: entry for entry in snapshot_entries}
    require(len(snapshot_index) == len(snapshot_entries), "Duplicate raw snapshot paths")
    require(snapshots["snapshots"] == snapshot_entries, "Published snapshot manifest diverges from raw manifest")
    validate_provenance(records, snapshot_index)
    for protein in proteins:
        validate_provenance([{"record_id": protein["protein_id"], "provenance": protein["provenance"]}], snapshot_index)
    for entry in snapshot_entries:
        path = ROOT / entry["path"]
        require(path.is_file(), f"Missing raw snapshot: {entry['path']}")
        require(sha256_file(path) == entry["sha256"], f"Raw snapshot hash mismatch: {entry['path']}")

    license_ids = {item["license"] for item in corpus["licenses"]["sources"]}
    require(license_ids == {"CC-BY-4.0", "CC-BY-SA-3.0", "CC0-1.0"}, "Source license set mismatch")
    require(corpus["licenses"]["aggregate_license"] is None, "Mixed-source bundle has blanket license")
    license_expectations = {
        "raw/licenses/uniprot-license.txt": "CC BY 4.0",
        "raw/licenses/rhea-license.txt": "CC BY 4.0",
        "raw/licenses/chembl-license-extract.txt": "Creative Commons Attribution-Share Alike 3.0",
        "raw/licenses/pdbe-license-extract.txt": "CC0",
    }
    for rel_path, phrase in license_expectations.items():
        require(phrase in (ROOT / rel_path).read_text(encoding="utf-8"), f"License evidence changed: {rel_path}")

    legal_code = {
        "CC-BY-4.0": ("https://creativecommons.org/licenses/by/4.0/legalcode", "Attribution 4.0 International"),
        "CC-BY-SA-3.0": ("https://creativecommons.org/licenses/by-sa/3.0/legalcode", "Attribution-ShareAlike 3.0 Unported"),
        "CC0-1.0": ("https://creativecommons.org/publicdomain/zero/1.0/legalcode", "CC0 1.0 Universal"),
    }
    for license_id, (_url, title) in legal_code.items():
        text_path = ROOT / "LICENSES" / f"{license_id}.txt"
        require(text_path.is_file() and title in text_path.read_text(encoding="utf-8"), f"Missing legal code: {text_path}")
    for source in corpus["licenses"]["sources"]:
        require(source["official_license_url"] == legal_code[source["license"]][0], f"License URL mismatch: {source['source']}")
        require((ROOT / source["license_text"]).is_file(), f"Missing license text: {source['license_text']}")
    record_licenses = {
        "experimental_structure_mapping": "CC-BY-4.0",
        "assay_activity": "CC-BY-SA-3.0",
        "reaction_assignment": "CC-BY-4.0",
    }
    for record in records:
        expected = record_licenses[record["evidence_type"]]
        require(record.get("license") == expected, f"Record license mismatch: {record['record_id']}")
        require(record.get("license_url") == legal_code[expected][0], f"Record license URL mismatch: {record['record_id']}")
    require(activities_file.get("license") == "CC-BY-SA-3.0", "activities.json lacks its CC BY-SA 3.0 license")
    require(activities_file.get("license_url") == legal_code["CC-BY-SA-3.0"][0], "activities.json license URL mismatch")
    require("10.6019/CHEMBL.database." in activities_file.get("citation", ""), "activities.json lacks the ChEMBL citation")

    validate_checksum_file(ROOT / "raw")
    validate_checksum_file(ROOT / "corpus")
    print(
        "OK: 15 proteins, 80 exact experimental structures, 28 direct assay activities, "
        "12 UniProt/Rhea assignments, 120 total records; provenance and checksums valid"
    )


if __name__ == "__main__":
    main()
