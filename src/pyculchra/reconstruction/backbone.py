"""Y7 backbone reconstruction extracted from Yitian_Exercises.ipynb.

Implementation: Yitian (AI-assisted); human review pending.
"""
import numpy as np
from ..geometry import fit_rigid_transform

def _signed_ca_distance(ca1, ca2, ca3, ca4):
    """Distance ca1-ca4, negative when the four CA atoms are left-handed."""
    turn = np.dot(np.cross(ca2 - ca1, ca3 - ca2), ca4 - ca3)
    distance = np.linalg.norm(ca4 - ca1)
    return -distance if turn < 0 else distance


def _extend_ca_trace(ca):
    """Add the two missing CA positions at each end (PULCHRA prepare_rbins).

    The source superimposes the first/second three CA atoms and moves the
    five-atom block, which is a rigid extrapolation of the chain direction.
    """
    ca = np.asarray(ca, dtype=float)
    if len(ca) < 5:
        raise ValueError('PULCHRA needs at least 5 CA atoms per chain')
    extended = np.empty((len(ca) + 4, 3))
    extended[2:-2] = ca
    rotation, shift = fit_rigid_transform(ca[2:5], ca[0:3])
    head = ca[0:5] @ rotation + shift
    extended[0], extended[1] = head[0], head[1]
    rotation, shift = fit_rigid_transform(ca[-5:-2], ca[-3:])
    tail = ca[-5:] @ rotation + shift
    extended[-2], extended[-1] = tail[3], tail[4]
    return extended


def _clamp(value, low, high):
    return max(low, min(high, value))


def _trace_bins(extended):
    """PULCHRA descriptor bins for every four-CA window (n + 1 windows)."""
    bins = []
    for i in range(len(extended) - 3):
        ca1, ca2, ca3, ca4 = extended[i:i + 4]
        r13_1 = np.linalg.norm(ca3 - ca1)
        r13_2 = np.linalg.norm(ca4 - ca2)
        r14 = _signed_ca_distance(ca1, ca2, ca3, ca4)
        bins.append((
            _clamp(int((r13_1 - 4.6) / 0.3), 0, 9),
            _clamp(int((r13_2 - 4.6) / 0.3), 0, 9),
            _clamp(int((r14 + 11.0) / 0.3), 0, 73),
        ))
    return bins


def _lookup_fragment(table, target_bins):
    """Best table fragment; the C loop stops at the sentinel or an exact hit."""
    best_index = None
    best_score = 1000.0
    for index, entry in enumerate(table):
        if entry['bins'][0] < 0:
            break
        score = (abs(entry['bins'][0] - target_bins[0])
                 + abs(entry['bins'][1] - target_bins[1])
                 + 0.2 * abs(entry['bins'][2] - target_bins[2]))
        if score < best_score:
            best_index, best_score = index, score
        if score <= 1e-3:
            break
    if best_index is None:
        raise ValueError('backbone library has no usable fragment')
    return best_index, best_score


def rebuild_backbone(ca_record, backbone_library, options):
    """Rebuild N/C/O (and keep CA) from a CA trace using PULCHRA 3.07 tables.

    The window CA[i-2..i+1] selects a fragment that supplies C and O of
    residue i-1 and N of residue i, exactly as rebuild_backbone() in
    pulchra.c. Atoms are never connected across a segment break.
    """
    keys = [tuple(key) for key in ca_record.get('atom_keys', [])]
    names = list(ca_record.get('residue_names', []))
    coords = np.asarray(ca_record.get('coords', np.empty((0, 3))), dtype=float)
    segments = [list(segment) for segment in ca_record.get('segments', [])]
    options = dict(options or {})

    unknown = set(options) - {'include_ca', 'library_version', 'source_commit'}
    if unknown:
        raise ValueError('unsupported options: %s' % sorted(unknown))
    if coords.shape != (len(keys), 3) or len(names) != len(keys):
        raise ValueError('atom_keys, residue_names and coords must agree')

    if not keys:
        return {'atom_keys': [], 'residue_names': [], 'coords': np.empty((0, 3)),
                'segments': [], 'diagnostics': []}
    if not backbone_library:
        raise ValueError('a reviewed backbone library is required')
    if not segments:
        raise ValueError('a nonempty CA record needs explicit segments')
    flat = [row for segment in segments for row in segment]
    if flat != list(range(len(keys))):
        raise ValueError('segments must list the CA rows once, in chain order')

    include_ca = bool(options.get('include_ca', True))
    atom_order = ['N', 'CA', 'C', 'O', 'OXT']
    rebuilt = {}
    diagnostics = []
    out_pieces = []
    out_segments = []

    for segment_number, segment in enumerate(segments):
        if len(segment) < 5:
            for row in segment:
                diagnostics.append({
                    'residue_key': keys[row],
                    'segment': segment_number,
                    'status': 'excluded',
                    'reason': 'PULCHRA needs at least 5 CA atoms per segment',
                })
            continue

        trace = coords[segment]
        extended = _extend_ca_trace(trace)
        bins = _trace_bins(extended)

        for window in range(len(bins)):
            is_proline = window > 0 and names[segment[window - 1]] == 'PRO'
            table_name = 'proline' if is_proline else 'standard'
            table = backbone_library[table_name]
            library_index, score = _lookup_fragment(table, bins[window])
            fragment = table[library_index]['coords']

            rotation, shift = fit_rigid_transform(fragment[:4], extended[window:window + 4])
            placed = fragment @ rotation + shift
            fitted_ca = fragment[:4] @ rotation + shift
            fit_rmsd = float(np.sqrt(np.mean(np.sum(
                (extended[window:window + 4] - fitted_ca) ** 2, axis=1))))

            if window > 0:
                previous = segment[window - 1]
                rebuilt[(previous, 'C')] = placed[4]
                rebuilt[(previous, 'O')] = placed[5]
            if window < len(segment):
                rebuilt[(segment[window], 'N')] = placed[6]
            else:
                rebuilt[(segment[window - 1], 'OXT')] = placed[6]

            diagnostics.append({
                'residue_key': keys[segment[min(window, len(segment) - 1)]],
                'segment': segment_number,
                'window': window,
                'bins': bins[window],
                'library': table_name,
                'library_index': library_index,
                'library_bins': tuple(table[library_index]['bins']),
                'matching_score': score,
                'fit_rmsd': fit_rmsd,
                'status': 'rebuilt',
                'reason': None,
            })

        start = len(out_pieces)
        for row in segment:
            for atom in atom_order:
                if atom == 'CA':
                    if include_ca:
                        out_pieces.append((row, atom))
                elif (row, atom) in rebuilt:
                    out_pieces.append((row, atom))
        out_segments.append(list(range(start, len(out_pieces))))

    out_keys, out_names, out_coords = [], [], []
    for row, atom in out_pieces:
        out_keys.append(keys[row][:3] + (atom,))
        out_names.append(names[row])
        out_coords.append(coords[row] if atom == 'CA' else rebuilt[(row, atom)])

    return {
        'atom_keys': out_keys,
        'residue_names': out_names,
        'coords': np.array(out_coords, dtype=float).reshape(-1, 3),
        'segments': out_segments,
        'diagnostics': diagnostics,
        'library_version': options.get('library_version', backbone_library.get('version')),
        'source_commit': options.get('source_commit'),
    }
