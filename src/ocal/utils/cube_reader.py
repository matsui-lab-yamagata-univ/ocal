"""
Cube file reader and orbital metrics calculator.

Computes the orbital tail fraction (OTF) and the participation ratio per
molecular volume (PR/V) from Gaussian cubegen MO cube files.

Only cubes with a Bohr-unit (all-positive voxel counts), orthogonal grid
are supported. Å-unit and non-orthogonal grids are rejected when the cube
is read, so neither OTF nor PR/V is computed for them.
"""
from __future__ import annotations

import csv
import math
from itertools import product
from pathlib import Path
from typing import NamedTuple

import numpy as np
from numpy import linalg as LA

from .vdw_volume import calc_vdw_union_volume

a0 = 0.5291772083

ELEMENT_PROP_PATH = (
    Path(__file__).resolve().parent.parent / "constants" / "element_properties.csv"
)

_ORTHO_TOL = 1e-6


def load_element_tables(
    path: Path = ELEMENT_PROP_PATH,
) -> tuple[dict[str, float], dict[str, int], dict[int, str]]:
    """
    Load vdW radii and atomic-number tables from element_properties.csv.

    Reads with ``utf-8-sig`` so a leading BOM does not corrupt the first
    column name. Rows with ``number == 0`` (Dummy ``X``) are omitted from
    the atomic-number maps. Empty, non-finite, or non-positive
    ``vdw_radius`` values are omitted from the radius map (lookups then
    fail explicitly). For ``number → symbol``, the first CSV row wins
    (so ``1 → "H"``, not ``"D"``).

    Parameters
    ----------
    path : Path
        Path to ``element_properties.csv``.

    Returns
    -------
    tuple[dict[str, float], dict[str, int], dict[int, str]]
        ``(van_der_waals_radius, elements_num, rev_elements_num)``.
    """
    van_der_waals_radius: dict[str, float] = {}
    elements_num: dict[str, int] = {}
    rev_elements_num: dict[int, str] = {}

    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                number = int(row["number"])
            except (KeyError, TypeError, ValueError):
                continue

            symbol = (row.get("symbol") or "").strip()
            if not symbol:
                continue

            if number != 0:
                elements_num[symbol] = number
                if number not in rev_elements_num:
                    rev_elements_num[number] = symbol

            try:
                radius = float(row["vdw_radius"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(radius) and radius > 0.0):
                continue
            van_der_waals_radius[symbol] = radius

    return van_der_waals_radius, elements_num, rev_elements_num


van_der_waals_radius, elements_num, rev_elements_num = load_element_tables()


class CubeData(NamedTuple):
    """Parsed contents of a cube file."""

    #: Voxel coordinates in Å, shape ``(n_voxels, 3)``, or ``(0, 3)``
    #: when coordinates were not requested.
    voxel_coords: np.ndarray
    #: Normalized probability density per voxel, shape ``(n_voxels,)``.
    norm_voxel_val: np.ndarray
    #: Atomic symbols, shape ``(n_atoms,)``.
    single_symbols: np.ndarray
    #: Atomic coordinates in Å, shape ``(n_atoms, 3)``.
    single_coords: np.ndarray
    #: vdW radii in Å, shape ``(n_atoms,)``.
    single_VDWrad: np.ndarray
    #: Voxel volume ``dV`` in Å³.
    voxel_volume_ang3: float


class CubeMetrics(NamedTuple):
    """OTF and PR/V results for one cube file."""

    #: Orbital tail fraction.
    otf: float
    #: Participation ratio divided by the vdW union volume (dimensionless).
    pr_over_volume: float
    #: Participation ratio in Å³.
    pr_ang3: float
    #: vdW union volume of the molecule in Å³.
    vdw_volume_ang3: float
    #: Sum of the normalized density over the whole grid (should be ~1).
    density_sum: float


def _parse_grid_header(
    header_lines: list[str], filepath: Path
) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse and validate the grid lines (lines 4–6) of a cube header.

    Parameters
    ----------
    header_lines : list[str]
        The three grid lines, each ``(n_i, v_ix, v_iy, v_iz)``.
    filepath : Path
        Cube file path, used in error messages.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        ``(voxel_num, voxel_interval)``: voxel counts of shape ``(3,)``
        and grid vectors in Bohr of shape ``(3, 3)``.

    Raises
    ------
    ValueError
        If a line cannot be parsed, voxel counts are zero or mixed in
        sign, the grid is an Å-unit (all-negative) cube, or the grid
        vectors are not orthogonal.
    """
    if len(header_lines) != 3:
        raise ValueError(f"{filepath}: truncated header (need lines 4–6)")

    counts: list[int] = []
    vectors: list[list[float]] = []
    for i, line in enumerate(header_lines, start=4):
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(
                f"{filepath}: line {i} must have 4 numeric fields, "
                f"got {len(parts)}"
            )
        try:
            count = int(parts[0])
            vec = [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError as exc:
            raise ValueError(
                f"{filepath}: line {i} is not parseable as "
                f"(count, vx, vy, vz): {line!r}"
            ) from exc
        counts.append(count)
        vectors.append(vec)

    if any(c == 0 for c in counts):
        raise ValueError(f"{filepath}: voxel count must be non-zero, got {counts}")
    if any(c > 0 for c in counts) and any(c < 0 for c in counts):
        raise ValueError(f"{filepath}: mixed-sign voxel counts {counts}")
    if all(c < 0 for c in counts):
        raise ValueError(
            f"{filepath}: Å-unit cube grids (all-negative voxel counts) "
            f"are not supported, got {counts}"
        )

    mat = np.asarray(vectors, dtype=float)
    for i in range(3):
        for j in range(i + 1, 3):
            dot = float(np.dot(mat[i], mat[j]))
            if abs(dot) > _ORTHO_TOL:
                raise ValueError(
                    f"{filepath}: non-orthogonal grid vectors "
                    f"(axes {i},{j} dot={dot})"
                )

    return np.asarray(counts, dtype=np.int64), mat


def read_cube(
    filepath: Path,
    plus_one: bool = True,
    with_coords: bool = True,
) -> CubeData:
    """
    Read a Gaussian MO cube file and validate its grid.

    Only Bohr-unit (all-positive voxel counts), orthogonal grids are
    accepted. All lengths in the returned data are in Å.

    Parameters
    ----------
    filepath : Path
        The path to the cube file.
    plus_one : bool
        If ``True`` (default), skip the MO-index line that follows the
        atom block (present when the atom count is negative, as in
        cubegen MO cubes).
    with_coords : bool
        If ``True`` (default), build the full ``(n_voxels, 3)`` voxel
        coordinate array. If ``False``, skip that construction and
        return an empty ``np.empty((0, 3))`` instead.

    Returns
    -------
    CubeData
        Parsed voxel data, atoms, and voxel volume.

    Raises
    ------
    ValueError
        If the header is invalid or unsupported (see
        ``_parse_grid_header``), an atomic number or vdW radius is
        missing, or the number of voxel values does not match the header
        dimensions.
    """
    with open(filepath, "r", encoding="utf-8") as reader:
        lines = reader.readlines()

    if len(lines) < 6:
        raise ValueError(f"{filepath}: truncated header (need lines 1–6)")

    voxel_num, voxel_interval_bohr = _parse_grid_header(lines[3:6], filepath)

    # Third line: atom count and origin. Bohr units are guaranteed by
    # _parse_grid_header (all-positive voxel counts), so the origin is
    # always converted; the sign of the atom count only flags the MO line.
    third_line = lines[2].split()
    atoms_num = abs(int(third_line[0]))
    voxel_origin = np.array([float(i) for i in third_line[1:4]]) * a0
    voxel_interval = voxel_interval_bohr * a0

    if with_coords:
        voxel_mut = np.array(
            list(product(range(voxel_num[0]), range(voxel_num[1]), range(voxel_num[2])))
        )
        voxel_coords = voxel_origin + (voxel_mut @ voxel_interval)
    else:
        voxel_coords = np.empty((0, 3))

    expected = int(voxel_num[0] * voxel_num[1] * voxel_num[2])
    voxel_val = np.empty((expected,))

    if plus_one:
        vol_start = 5 + atoms_num + 2
    else:
        vol_start = 5 + atoms_num + 1

    single_symbols = np.empty((atoms_num,), dtype="object")
    single_coords = np.empty((atoms_num, 3))
    single_VDWrad = np.empty((atoms_num,))
    for ind, line in enumerate(lines[6 : 6 + atoms_num]):
        line_split = line.strip().split()
        symbol_ind = int(line_split[0])
        try:
            symbol = rev_elements_num[symbol_ind]
        except KeyError as exc:
            raise ValueError(f"{filepath}: unknown atomic number {symbol_ind}") from exc
        try:
            radius = van_der_waals_radius[symbol]
        except KeyError as exc:
            raise ValueError(
                f"{filepath}: no vdw_radius in element_properties.csv for {symbol!r}"
            ) from exc
        single_symbols[ind] = symbol
        single_VDWrad[ind] = radius
        single_coords[ind] = np.array([float(i) for i in line_split[2:]]) * a0

    count = 0
    for line in lines[vol_start:]:
        line_split = [float(i) for i in line.strip().split()]
        n = len(line_split)
        if n == 0:
            continue
        end = count + n
        if end <= expected:
            voxel_val[count:end] = np.array(line_split)
        elif count < expected:
            voxel_val[count:expected] = np.array(line_split[: expected - count])
        count = end

    if count != expected:
        raise ValueError(
            f"{filepath}: expected {expected} voxel values, got {count}"
        )

    # Orthogonal grid: product of vector norms equals |det(V)|.
    dV_bohr3 = float(
        LA.norm(voxel_interval_bohr[0])
        * LA.norm(voxel_interval_bohr[1])
        * LA.norm(voxel_interval_bohr[2])
    )
    norm_voxel_val = voxel_val**2 * dV_bohr3

    return CubeData(
        voxel_coords=voxel_coords,
        norm_voxel_val=norm_voxel_val,
        single_symbols=single_symbols,
        single_coords=single_coords,
        single_VDWrad=single_VDWrad,
        voxel_volume_ang3=abs(float(np.linalg.det(voxel_interval_bohr))) * a0**3,
    )


def read_cubefile(
    filepath: Path,
    plus_one: bool = True,
    with_coords: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a cube file and return voxel and atom arrays.

    Thin wrapper around ``read_cube`` kept for backward compatibility.

    Parameters
    ----------
    filepath : Path
        The path to the cube file.
    plus_one : bool
        If ``True`` (default), skip the MO-index line after the atom block.
    with_coords : bool
        If ``True`` (default), build the voxel coordinate array.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        ``(voxel_coords, norm_voxel_val, single_symbols, single_coords,
        single_VDWrad)``.

    Raises
    ------
    ValueError
        See ``read_cube``.
    """
    data = read_cube(filepath, plus_one=plus_one, with_coords=with_coords)
    return (
        data.voxel_coords,
        data.norm_voxel_val,
        data.single_symbols,
        data.single_coords,
        data.single_VDWrad,
    )


def read_voxel_volume(cube_path: Path) -> float:
    """
    Parse the cube header grid vectors and return the voxel volume in Å³.

    Parameters
    ----------
    cube_path : Path
        Path to a Gaussian MO cube file.

    Returns
    -------
    float
        Voxel volume ``dV`` in Å³.

    Raises
    ------
    ValueError
        If the header is truncated, invalid, or unsupported (see
        ``_parse_grid_header``).
    """
    with cube_path.open(encoding="utf-8") as fh:
        for _ in range(3):
            if not fh.readline():
                raise ValueError(f"{cube_path}: truncated header (need lines 4–6)")
        header_lines = [line for line in (fh.readline() for _ in range(3)) if line]

    _counts, mat = _parse_grid_header(header_lines, cube_path)
    return abs(float(np.linalg.det(mat))) * a0**3


def calc_distance_matrix(
    voxel_coords: np.ndarray,
    single_coords: np.ndarray,
    single_VDWrad: np.ndarray,
) -> np.ndarray:
    """
    Calculate the signed distance from each voxel to the vdW surface.

    Parameters
    ----------
    voxel_coords : np.ndarray
        The voxel coordinates.
    single_coords : np.ndarray
        The coordinates of single atoms.
    single_VDWrad : np.ndarray
        The van der Waals radii of single atoms.

    Returns
    -------
    np.ndarray
        Minimum over atoms of ``|r_voxel - r_atom| - r_vdW``.

    Raises
    ------
    ValueError
        If ``single_coords`` contains no atoms.
    """
    distance_matrix: np.ndarray | None = None
    for ind in range(single_coords.shape[0]):
        temp_distance_matrix = (
            LA.norm(voxel_coords - single_coords[ind], axis=-1).reshape(1, -1)
            - single_VDWrad[ind]
        )
        if ind == 0:
            distance_matrix = temp_distance_matrix.copy()
        else:
            distance_matrix = np.min(
                np.concatenate([distance_matrix, temp_distance_matrix], axis=0),
                axis=0,
                keepdims=True,
            )
    if distance_matrix is None:
        raise ValueError("single_coords must contain at least one atom")
    return distance_matrix.reshape(-1)


def _otf_from_data(data: CubeData) -> float:
    """
    Compute OTF from parsed cube data.

    Parameters
    ----------
    data : CubeData
        Parsed cube with voxel coordinates.

    Returns
    -------
    float
        Sum of normalized density over voxels outside the vdW surface.
    """
    distance_matrix = calc_distance_matrix(
        data.voxel_coords, data.single_coords, data.single_VDWrad
    )
    return float(np.sum(data.norm_voxel_val[distance_matrix >= 0]))


def _pr_from_data(
    data: CubeData,
    cube_path: Path,
    spacing: float,
    vdw_volume_ang3: float | None,
) -> tuple[float, float, float]:
    """
    Compute PR/V, PR, and vdW volume from parsed cube data.

    Uses the normalization-invariant form
    ``pr_ang3 = dV_ang3 * (Σ p_i)² / Σ p_i²``.

    Parameters
    ----------
    data : CubeData
        Parsed cube data.
    cube_path : Path
        Cube file path, used in error messages.
    spacing : float
        Grid spacing in Å for the vdW union volume.
    vdw_volume_ang3 : float | None
        Optional precomputed vdW union volume in Å³.

    Returns
    -------
    tuple[float, float, float]
        ``(pr_over_volume, pr_ang3, vdw_volume_ang3)``.

    Raises
    ------
    ValueError
        If ``Σ p_i²`` is non-positive, a provided volume is invalid, or a
        result is non-finite / non-positive.
    """
    density_sum = float(np.sum(data.norm_voxel_val))
    denom = float(np.sum(data.norm_voxel_val**2))
    if not (denom > 0.0):
        raise ValueError(f"{cube_path}: Σ p_i² must be > 0, got {denom}")
    pr_ang3 = data.voxel_volume_ang3 * density_sum**2 / denom

    if vdw_volume_ang3 is not None:
        if not (math.isfinite(vdw_volume_ang3) and vdw_volume_ang3 > 0.0):
            raise ValueError(
                f"{cube_path}: vdw_volume_ang3 must be finite and > 0, "
                f"got {vdw_volume_ang3}"
            )
        volume = float(vdw_volume_ang3)
    else:
        volume = calc_vdw_union_volume(
            data.single_coords, data.single_VDWrad, spacing=spacing
        )

    pr_over_volume = pr_ang3 / volume

    if not (math.isfinite(pr_ang3) and pr_ang3 > 0.0):
        raise ValueError(
            f"{cube_path}: pr_ang3 must be finite and > 0, got {pr_ang3}"
        )
    if not (math.isfinite(pr_over_volume) and pr_over_volume > 0.0):
        raise ValueError(
            f"{cube_path}: pr_over_volume must be finite and > 0, "
            f"got {pr_over_volume}"
        )
    return pr_over_volume, pr_ang3, volume


def _check_density_sum(density_sum: float, cube_path: Path) -> None:
    """
    Validate that the density sum is finite.

    Parameters
    ----------
    density_sum : float
        Sum of the normalized density.
    cube_path : Path
        Cube file path, used in error messages.

    Raises
    ------
    ValueError
        If ``density_sum`` is not finite.
    """
    if not math.isfinite(density_sum):
        raise ValueError(
            f"{cube_path}: density_sum must be finite, got {density_sum}"
        )


def calc_otf(cube_path: Path) -> tuple[float, float]:
    """
    Compute the orbital tail fraction and density sum for a cube file.

    The fraction is the sum of normalized probability density over voxels
    outside the van der Waals surface (signed distance ``>= 0``).

    Parameters
    ----------
    cube_path : Path
        Path to a Gaussian MO cube file (Bohr-unit, orthogonal grid).

    Returns
    -------
    tuple[float, float]
        ``(otf, density_sum)``.

    Raises
    ------
    ValueError
        If the cube is invalid or unsupported (see ``read_cube``).
    """
    data = read_cube(cube_path, with_coords=True)
    density_sum = float(np.sum(data.norm_voxel_val))
    _check_density_sum(density_sum, cube_path)
    return _otf_from_data(data), density_sum


def calc_pr(
    cube_path: Path,
    spacing: float = 0.2,
    vdw_volume_ang3: float | None = None,
) -> tuple[float, float, float, float]:
    """
    Compute normalized participation ratio and related quantities.

    Parameters
    ----------
    cube_path : Path
        Path to a Gaussian MO cube file (Bohr-unit, orthogonal grid).
    spacing : float
        Grid spacing in Å for the vdW union volume (default 0.2).
    vdw_volume_ang3 : float | None
        Optional precomputed molecular vdW union volume in Å³. When
        given (finite and ``> 0``), it is reused instead of recomputed.

    Returns
    -------
    tuple[float, float, float, float]
        ``(pr_over_volume, pr_ang3, vdw_volume_ang3, density_sum)``.

    Raises
    ------
    ValueError
        If the cube is invalid or unsupported, or a result is invalid.
    """
    data = read_cube(cube_path, with_coords=False)
    density_sum = float(np.sum(data.norm_voxel_val))
    _check_density_sum(density_sum, cube_path)
    pr_over_volume, pr_ang3, volume = _pr_from_data(
        data, cube_path, spacing, vdw_volume_ang3
    )
    return pr_over_volume, pr_ang3, volume, density_sum


def calc_otf_pr(
    cube_path: Path,
    spacing: float = 0.2,
    vdw_volume_ang3: float | None = None,
) -> CubeMetrics:
    """
    Compute OTF and PR/V for a cube file, reading it only once.

    The cube grid is validated before either metric is computed, so an
    unsupported cube yields neither value.

    Parameters
    ----------
    cube_path : Path
        Path to a Gaussian MO cube file (Bohr-unit, orthogonal grid).
    spacing : float
        Grid spacing in Å for the vdW union volume (default 0.2).
    vdw_volume_ang3 : float | None
        Optional precomputed molecular vdW union volume in Å³.

    Returns
    -------
    CubeMetrics
        OTF, PR/V, PR, vdW volume, and density sum.

    Raises
    ------
    ValueError
        If the cube is invalid or unsupported, or a result is invalid.
    """
    data = read_cube(cube_path, with_coords=True)
    density_sum = float(np.sum(data.norm_voxel_val))
    _check_density_sum(density_sum, cube_path)
    otf = _otf_from_data(data)
    pr_over_volume, pr_ang3, volume = _pr_from_data(
        data, cube_path, spacing, vdw_volume_ang3
    )
    return CubeMetrics(
        otf=otf,
        pr_over_volume=pr_over_volume,
        pr_ang3=pr_ang3,
        vdw_volume_ang3=volume,
        density_sum=density_sum,
    )
