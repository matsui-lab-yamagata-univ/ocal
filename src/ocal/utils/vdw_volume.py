"""Van der Waals union volume via cubic-grid integration."""
from __future__ import annotations

import math

import numpy as np

_MAX_GRID_POINTS = 200_000_000


def calc_vdw_union_volume(
    coords: np.ndarray,
    radii: np.ndarray,
    spacing: float = 0.2,
) -> float:
    """
    Approximate the union volume of hard vdW spheres on a cubic grid.

    Each atom marks only the grid points inside its own bounding box,
    so the cost scales with Σ(sphere volume) / spacing³ rather than
    n_atoms × n_grid.

    Parameters
    ----------
    coords : np.ndarray
        Atom coordinates in Å, shape ``(n, 3)``.
    radii : np.ndarray
        vdW radii in Å, shape ``(n,)``.
    spacing : float
        Isotropic grid spacing in Å (default 0.2).

    Returns
    -------
    float
        Approximate union volume in Å³.

    Raises
    ------
    ValueError
        If inputs are empty/invalid, the grid exceeds the size limit,
        or the computed volume is not finite and positive.
    """
    coords = np.asarray(coords, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(
            f"coords must have shape (n, 3), got {coords.shape}"
        )
    if radii.ndim != 1 or radii.shape[0] != coords.shape[0]:
        raise ValueError(
            f"radii must have shape ({coords.shape[0]},), got {radii.shape}"
        )
    if coords.shape[0] == 0:
        raise ValueError("no atoms to compute volume for")
    if not (math.isfinite(spacing) and spacing > 0.0):
        raise ValueError(f"spacing must be > 0, got {spacing}")
    if not np.all(np.isfinite(coords)) or not np.all(np.isfinite(radii)):
        raise ValueError("coords and radii must be finite")
    if not np.all(radii > 0.0):
        raise ValueError("all radii must be > 0")

    lo = np.min(coords - radii[:, None], axis=0)
    hi = np.max(coords + radii[:, None], axis=0)
    # Inclusive grid covering [lo, hi]
    n_xyz = np.floor((hi - lo) / spacing).astype(int) + 1
    n_points = int(np.prod(n_xyz))
    if n_points > _MAX_GRID_POINTS:
        raise ValueError(
            f"grid too large: {n_points} points "
            f"(limit {_MAX_GRID_POINTS}); bbox={lo}..{hi}, "
            f"spacing={spacing}"
        )
    if n_points <= 0:
        raise ValueError(f"degenerate grid shape {tuple(n_xyz)}")

    mark = np.zeros(tuple(n_xyz.tolist()), dtype=bool)
    for i in range(coords.shape[0]):
        c = coords[i]
        r = float(radii[i])
        r2 = r * r
        i_lo = np.maximum(
            np.ceil((c - r - lo) / spacing).astype(int), 0
        )
        i_hi = np.minimum(
            np.floor((c + r - lo) / spacing).astype(int),
            n_xyz - 1,
        )
        if np.any(i_lo > i_hi):
            continue
        xs = lo[0] + np.arange(i_lo[0], i_hi[0] + 1) * spacing
        ys = lo[1] + np.arange(i_lo[1], i_hi[1] + 1) * spacing
        zs = lo[2] + np.arange(i_lo[2], i_hi[2] + 1) * spacing
        dx = xs - c[0]
        dy = ys - c[1]
        dz = zs - c[2]
        d2 = (
            dx[:, None, None] ** 2
            + dy[None, :, None] ** 2
            + dz[None, None, :] ** 2
        )
        mark[
            i_lo[0] : i_hi[0] + 1,
            i_lo[1] : i_hi[1] + 1,
            i_lo[2] : i_hi[2] + 1,
        ] |= d2 <= r2

    volume = float(mark.sum()) * spacing**3
    if not (math.isfinite(volume) and volume > 0.0):
        raise ValueError(
            f"vdw_volume_ang3 must be finite and > 0, got {volume}"
        )
    return volume
