# Security

Report vulnerabilities privately through [GitHub private vulnerability reporting](https://github.com/draconitic-labs/carbonic-anhydrase-corpus/security/advisories/new), or at [support@draconitic.com](mailto:support@draconitic.com).

Fixes go to the `main` branch. There are no maintained older branches.

`make build` and `make offline` read only the snapshots under `raw/` and make no network requests when every snapshot is present. `make refresh` downloads new snapshots over HTTPS from the official UniProt, Rhea, ChEMBL and PDBe APIs. Every snapshot's SHA-256 is recorded in `raw/checksums.sha256` and checked by `make validate`.
