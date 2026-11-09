# Human alpha-carbonic anhydrase evidence corpus

This is a release-pinned evidence corpus for the 15 core
human alpha-carbonic anhydrase paralogs. It contains exactly 120 normalized
evidence records:

- 80 experimental PDB/SIFTS structure mappings;
- 28 ChEMBL assay/activity records with direct, exact single-protein targets;
- 12 UniProt/Rhea carbonic-anhydrase reaction assignments.

Raw official-API snapshots are kept with their URLs, retrieval metadata,
release metadata and SHA-256 checksums.

## Licences

- **Code** (`scripts/`, `tests/`, `schema/`, `Makefile`): Apache-2.0, see `LICENSE`.
- **Data** (`corpus/`, `raw/`): each evidence record carries `license` and `license_url`.
  - CC BY-SA 3.0: `corpus/activities.json`, the `assay_activity` records in
    `corpus/corpus.json`, and `raw/chembl/**` (ChEMBL release 37,
    doi:10.6019/CHEMBL.database.37).
  - CC BY 4.0: protein dossiers, structure mappings and reaction assignments,
    which copy UniProtKB entries, cross-references and annotations, and the
    UniProt and Rhea snapshots.
  - CC0 1.0: PDBe and wwPDB data under `raw/pdbe/`.

The full legal codes are in `LICENSES/`; `NOTICE` and `corpus/licenses.json`
give the attribution for each source.

## Build and verify

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
make build
make test
```

`make offline` rebuilds from the checked-in source snapshots without network
access. `make refresh` deliberately replaces those snapshots from the live
official APIs and may change the resulting release. A normal build reuses any
existing snapshot and only downloads a missing one.

The primary integration artifact is `corpus/corpus.json`. The same evidence
records are split into `structures.json`, `activities.json`, and
`reactions.json` for simpler consumers. `family.json` holds the 15 normalized
protein dossiers, and `source-snapshots.json` plus `licenses.json` describe
provenance and reuse terms.

## Selection policy

Structure quotas prevent CA2's very large structure literature from erasing
the rest of the family. Within each protein, UniProt PDB cross-references are
bucketed by experimental method and reported resolution band, sorted by PDB
identifier, then selected round-robin across buckets. Every selected entry is
required to map back to the exact UniProt accession through PDBe SIFTS.

ChEMBL records are selected independently for each mapped family target. The
selector starts from confidence-9 direct assays whose target is an exact human
`SINGLE PROTEIN` component for the pinned UniProt accession, ordered by assay
identifier. For each target it reads the first 50 such assays, inspects the
first (quota + 10) that pass the direct-target check, and takes the earliest
admissible activity identifier within each inspected assay. Candidates must have a reported standardized value, unit, and relation;
must not be flagged as a potential duplicate or carry a data-validity comment;
and must resolve to an exact compound InChIKey. Compound forms retain their
exact ChEMBL molecule identifier and parent-form identifier. No value or
potency is used to rank the slice.

Reaction assignments are one per catalytically annotated paralog. Each record
requires the UniProt catalytic-activity comment to cite `RHEA:10748`, and the
equation is independently snapshotted from Rhea.

## Scientific boundary

CA8, CA10, and CA11 are included as acatalytic alpha-carbonic-anhydrase-related
proteins because UniProt explicitly states that they lack catalytic activity.
CA11 nevertheless has ChEMBL measurements assigned to a direct protein target;
those source records are retained but marked with
`CATALYTIC_STATUS_MISMATCH`, and are not evidence that CA11 catalyzes the Rhea
reaction.

Four core paralogs currently have no experimental PDB cross-reference in the
pinned UniProt release. Their zero-structure state remains visible in the
family dossiers. Structure counts are therefore a curated slice, not a claim
of equal structural coverage.

Assay descriptions and structured parameters are copied as source fields;
missing pH, temperature, construct, or protocol details remain `null` or
absent. Cross-assay comparison is not prequalified. A consumer should refuse a
comparison unless target, construct, parameter type, relation, unit, assay
format, conditions, and molecular form are compatible.

## Licensing

See [Licences](#licences) above, `NOTICE`, `LICENSES/` and `corpus/licenses.json`.
