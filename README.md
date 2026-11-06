# PyCULCHRA_eval

Course project for Python reconstruction of PULCHRA stages and evaluation.
Currently implemented here: Yitian Y4 restricted PDB parsing, Y5 rigid fitting,
Y6 local frames, and Y7 backbone reconstruction. Side-chain reconstruction
and the full pipeline are not yet included in this repository.

## Install and test

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

```python
from pyculchra.io import parse_pdb_atoms
from pyculchra.geometry import fit_rigid_transform, make_local_frame
```

All distances use Angstroms. Atom identity is
`(chain, residue_number, insertion_code, atom_name)`.
Rigid transforms use row vectors: `moving @ R + t`.
Local frames store axes as columns: `local @ axes.T + origin`.
PDB parsing returns a list of `{"key": ..., "xyz": ...}` records, not a full
pipeline structure record; an explicit adapter will be needed.

See [repository plan](docs/REPOSITORY_PLAN.md) for the proposed four-person
integration layout and [provenance](docs/PROVENANCE.md) for authorship and limits.

## Y7 backbone reconstruction

```python
from pyculchra.libraries import parse_backbone_library
from pyculchra.reconstruction import rebuild_backbone

# Run from repository root. The table is repository data, not wheel package data.
library = parse_backbone_library("third_party/pulchra/nco_data.h")
# ca_record: atom_keys, residue_names, coords (N,3), explicit segments
# result = rebuild_backbone(ca_record, library, {})
```

Each continuous segment needs at least five CA atoms; shorter segments are
excluded with diagnostics. Output includes N/CA/C/O and terminal OXT. Options:
include_ca (default True), library_version, source_commit. No native non-CA
coordinates enter reconstruction. Segment boundaries are supplied by the caller.

The tests include one recorded 1UBQ C reference in the original frame. This is
limited evidence, not whole-program equivalence. Second-protein and C-derived
terminal/gap validation remain pending. See data/fixtures/README.md and
third_party/pulchra/SOURCE.md for fixture and upstream provenance.
