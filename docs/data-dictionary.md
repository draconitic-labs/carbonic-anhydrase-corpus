# Data dictionary

## Stable identities

- Protein: UniProt primary accession + NCBI taxonomy identifier + sequence
  version and CRC64.
- Structure: PDB identifier + exact UniProt accession. Every PDBe/SIFTS entity
  and chain mapping is retained inside the record; consumers that operate at
  chain resolution must use `entity_id`, `chain_id`, and `struct_asym_id`.
- ChEMBL activity: ChEMBL release + integer `activity_id`.
- Compound: exact ChEMBL molecule identifier and standard InChIKey. Parent and
  salt forms are not silently collapsed.
- Reaction assignment: UniProt accession + Rhea master reaction identifier.

## Evidence record types

`experimental_structure_mapping` records contain the UniProt PDB
cross-reference, the matching SIFTS ranges and identity/coverage, PDBe entry
metadata, experimental metadata, and entity-level construct context. `null`
mutation flags or expression hosts mean the API did not provide a value.

`assay_activity` records contain target-component identity, assay context,
reported and standardized measurements, exact compound form, and document
metadata. Relation operators such as `<`, `>`, `<=`, and `>=` are preserved;
they are never rewritten as equality. `pChEMBL` is retained only when supplied
by ChEMBL and is not calculated locally.

`reaction_assignment` records contain the exact UniProt catalytic-activity
comment, evidence references, and the matching Rhea equation and EC class.

Every evidence record has `license` (an SPDX identifier) and `license_url` (the
Creative Commons legal code): `CC-BY-4.0` for structure mappings and reaction
assignments, `CC-BY-SA-3.0` for assay activities.

## Comparison refusal vocabulary

The corpus reserves these machine-readable codes for consumers:

- `TARGET_MISMATCH`
- `TARGET_SCOPE_NOT_SINGLE_PROTEIN`
- `TAXON_MISMATCH`
- `VARIANT_OR_CONSTRUCT_MISMATCH`
- `CATALYTIC_STATUS_MISMATCH`
- `PARAMETER_TYPE_MISMATCH`
- `CENSORED_RELATION`
- `UNSUPPORTED_UNIT`
- `ASSAY_FORMAT_MISMATCH`
- `CONDITION_MISSING`
- `CONDITION_MISMATCH`
- `MOLECULE_FORM_AMBIGUOUS`
- `STEREOCHEMISTRY_AMBIGUOUS`
- `DUPLICATE_DEPENDENCE`
- `REACTION_DIRECTION_MISMATCH`
- `PREDICTED_VS_EXPERIMENTAL`
- `PROVENANCE_INCOMPLETE`

The v1 corpus marks only directly evidenced, record-local conditions. It does
not attempt to infer protocol compatibility from free text.

