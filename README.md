# PyCULCHRA_eval

Course project for Python reconstruction of PULCHRA stages and evaluation.
Currently implemented here: Yitian Y4 restricted PDB parsing, Y5 rigid fitting,
and Y6 local frames. Backbone/side-chain reconstruction and the full pipeline
are not yet included in this repository.

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
