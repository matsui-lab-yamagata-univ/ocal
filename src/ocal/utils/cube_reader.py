"""Cube file reader and orbital tail fraction (OTF) calculator."""
from __future__ import annotations

import csv
import math
from itertools import product
from pathlib import Path

import numpy as np
from numpy import linalg as LA

a0 = 0.5291772083

ELEMENT_PROP_PATH = (
    Path(__file__).resolve().parent.parent / "constants" / "element_properties.csv"
)


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


def read_cubefile(
    filepath: Path, plus_one: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read and process a cube file.

    Parameters
    ----------
    filepath : Path
        The path of the cube file to be processed.
    plus_one : bool
        A flag to adjust the starting line for volume data.

    Returns
    -------
    tuple
        A tuple containing voxel coordinates, normalized voxel values, atomic symbols,
        atomic coordinates, and van der Waals radii.
    """
    with open(filepath, "r", encoding="utf-8") as reader:
        lines = reader.readlines()

    # Extract the number of atoms and voxel origin from the third line
    third_line = lines[2].split()
    atoms_num = int(third_line[0])
    voxel_origin = np.array([float(i) for i in third_line[1:4]])
    if atoms_num < 0:
        voxel_origin *= a0
        atoms_num *= -1

    # Extract voxel numbers and intervals
    voxel_num = np.array(
        [int(lines[3].split()[0]), int(lines[4].split()[0]), int(lines[5].split()[0])],
        dtype=np.int32,
    )
    voxel_interval = np.array(
        [
            [float(i) for i in lines[3].strip().split()[1:]],
            [float(i) for i in lines[4].strip().split()[1:]],
            [float(i) for i in lines[5].strip().split()[1:]],
        ]
    )
    voxel_interval_ori = voxel_interval.copy()
    voxel_interval[voxel_num >= 0] *= a0
    voxel_num = np.abs(voxel_num)

    # Calculate voxel coordinates
    voxel_mut = np.array(
        list(product(range(voxel_num[0]), range(voxel_num[1]), range(voxel_num[2])))
    )
    voxel_coords = voxel_origin + (voxel_mut @ voxel_interval)
    voxel_val = np.empty((voxel_num[0] * voxel_num[1] * voxel_num[2],))

    # Parse atom information
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

    # Parse voxel values
    count = 0
    for line in lines[vol_start:]:
        line_split = line.strip().split()
        line_split = [float(i) for i in line_split]
        voxel_val[count : count + len(line_split)] = np.array(line_split)
        count += len(line_split)

    # Normalize voxel values
    norm_voxel_val = voxel_val**2 * (
        LA.norm(voxel_interval_ori[0])
        * LA.norm(voxel_interval_ori[1])
        * LA.norm(voxel_interval_ori[2])
    )

    return voxel_coords, norm_voxel_val, single_symbols, single_coords, single_VDWrad


def calc_distance_matrix(
    voxel_coords: np.ndarray,
    single_coords: np.ndarray,
    single_VDWrad: np.ndarray,
) -> np.ndarray:
    """
    Calculate the distance matrix.

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
        The calculated distance matrix.
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


def calc_otf(cube_path: Path) -> tuple[float, float]:
    """
    Compute the orbital tail fraction and density sum for a cube file.

    The fraction is the sum of normalized probability density over voxels
    outside the van der Waals surface (signed distance ``>= 0``).
    ``density_sum`` is the total of all normalized voxel values and is
    used as a sanity check that the cube is properly normalized.

    Parameters
    ----------
    cube_path : Path
        Path to a Gaussian MO cube file.

    Returns
    -------
    tuple[float, float]
        ``(otf, density_sum)`` where ``otf`` is the orbital tail fraction
        and ``density_sum`` is ``norm_voxel_val.sum()``.
    """
    voxel_coords, norm_voxel_val, _symbols, single_coords, single_VDWrad = (
        read_cubefile(cube_path)
    )
    distance_matrix = calc_distance_matrix(
        voxel_coords, single_coords, single_VDWrad
    )
    otf = float(np.sum(norm_voxel_val[distance_matrix >= 0]))
    density_sum = float(np.sum(norm_voxel_val))
    return otf, density_sum
