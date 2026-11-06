"""Extracted from Yitian_Exercises.ipynb; see docs/PROVENANCE.md."""

import numpy as np


def parse_pdb_atoms(lines):
    """Parse the restricted Y4 ATOM format, preserving accepted atom order."""
    records = []
    seen = set()
    for number, line in enumerate(lines, start=1):
        if not isinstance(line, str):
            raise ValueError(f'Line {number}: expected text')
        record_type = line[:6].strip()
        if record_type == 'ENDMDL':
            break
        if record_type != 'ATOM':
            continue
        if len(line.rstrip('\r\n')) < 17:
            raise ValueError(f'Line {number}: truncated ATOM record')
        if line[16] not in (' ', 'A'):
            continue
        if len(line.rstrip('\r\n')) < 54:
            raise ValueError(f'Line {number}: missing fixed-width fields')
        atom_name = line[12:16].strip()
        if not atom_name:
            raise ValueError(f'Line {number}: missing atom name')
        try:
            residue_number = int(line[22:26])
            xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        except ValueError as exc:
            raise ValueError(f'Line {number}: invalid number field') from exc
        if not np.isfinite(xyz).all():
            raise ValueError(f'Line {number}: nonfinite coordinate')
        key = (line[21].strip(), residue_number, line[26].strip(), atom_name)
        if key in seen:
            raise ValueError(f'Line {number}: duplicate accepted atom {key}')
        seen.add(key)
        records.append({'key': key, 'xyz': xyz})
    return records
