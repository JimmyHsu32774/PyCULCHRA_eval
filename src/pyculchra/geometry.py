"""Extracted from Yitian_Exercises.ipynb; see docs/PROVENANCE.md."""

import numpy as np


def fit_rigid_transform(moving, reference):
    """Least-squares rotation R and translation t so that moving @ R + t matches reference."""
    moving = np.asarray(moving, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if moving.ndim != 2 or moving.shape != reference.shape or moving.shape[1] != 3:
        raise ValueError('need two matching (N, 3) arrays')
    if len(moving) < 3:
        raise ValueError('need at least three points')
    if not np.isfinite(moving).all() or not np.isfinite(reference).all():
        raise ValueError('coordinates must be finite')

    moving_mean = moving.mean(axis=0)
    reference_mean = reference.mean(axis=0)
    centred_moving = moving - moving_mean
    centred_reference = reference - reference_mean
    if (np.linalg.matrix_rank(centred_moving) < 2
            or np.linalg.matrix_rank(centred_reference) < 2):
        raise ValueError('points must span at least two dimensions')

    left, _, right = np.linalg.svd(centred_moving.T @ centred_reference)
    # a negative determinant would turn the data inside out; multiply one axis by -1
    flip = np.diag([1.0, 1.0, np.sign(np.linalg.det(left @ right))])
    rotation = left @ flip @ right
    translation = reference_mean - moving_mean @ rotation
    return rotation, translation


def make_local_frame(origin, forward, plane_point):
    """Right-handed frame; global coordinates are local @ axes.T + origin."""
    origin = np.asarray(origin, dtype=float)
    forward = np.asarray(forward, dtype=float)
    plane_point = np.asarray(plane_point, dtype=float)
    for point in (origin, forward, plane_point):
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError('expected three finite 3-vectors')

    e1 = forward - origin
    length = np.linalg.norm(e1)
    if length < 1e-12:
        raise ValueError('forward point must differ from origin')
    e1 = e1 / length

    e2 = plane_point - origin
    e2 = e2 - np.dot(e2, e1) * e1
    length = np.linalg.norm(e2)
    if length < 1e-12:
        raise ValueError('the three points must not be collinear')
    e2 = e2 / length

    e3 = np.cross(e1, e2)
    return origin.copy(), np.column_stack((e1, e2, e3))
