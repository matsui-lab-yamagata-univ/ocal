"""Structure file readers returning Cartesian coordinates in Angstrom."""
from __future__ import annotations

import functools
import re
from pathlib import Path

import numpy as np

from .cif_reader import CifReader
from .cube_reader import load_element_tables

print = functools.partial(print, flush=True)

_CHARGE_SPIN_RE = re.compile(r"^-?\d+\s+\d+$")
_ELEMENT_SYMBOL_RE = re.compile(r"^[A-Z][a-z]?$")

_van_der_waals_radius, _elements_num, rev_elements_num = load_element_tables()

CoordTuple = tuple[float, float, float]
Structure = tuple[list[str], list[CoordTuple]]


def read_gjf(path: Path) -> Structure:
    """
    Read atomic symbols and Cartesian coordinates from a Gaussian gjf/com file.

    Only the ``symbol x y z`` Cartesian format is supported. Charge and spin
    on the charge-spin line are ignored (always treated as ``0 1`` downstream).

    Parameters
    ----------
    path : Path
        Path to a ``.gjf`` or ``.com`` file.

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.

    Raises
    ------
    ValueError
        If no charge-spin line is found, no atoms are present, or a line uses
        an unsupported format (atomic number, freeze flag, Z-matrix, ONIOM).
    """
    symbols: list[str] = []
    coordinates: list[CoordTuple] = []

    with path.open(encoding="utf-8") as fh:
        lines = fh.readlines()

    start: int | None = None
    for i, line in enumerate(lines):
        if _CHARGE_SPIN_RE.match(line.strip()):
            start = i + 1
            break

    if start is None:
        raise ValueError(
            f"{path}: charge-spin line (e.g. '0 1') not found; "
            "only Cartesian 'symbol x y z' after charge-spin is supported"
        )

    for line in lines[start:]:
        stripped = line.strip()
        if not stripped:
            break
        parts = stripped.split()
        if len(parts) != 4:
            raise ValueError(
                f"{path}: unsupported atom line {stripped!r}; "
                "only 'symbol x y z' Cartesian format is supported "
                "(no atomic numbers, freeze flags, Z-matrix, or ONIOM)"
            )
        token, xs, ys, zs = parts
        if not _ELEMENT_SYMBOL_RE.match(token):
            raise ValueError(
                f"{path}: unsupported atom token {token!r}; "
                "only element-symbol Cartesian format is supported "
                "(no atomic numbers, freeze flags, Z-matrix, or ONIOM)"
            )
        try:
            coordinates.append((float(xs), float(ys), float(zs)))
        except ValueError as exc:
            raise ValueError(
                f"{path}: invalid coordinates in line {stripped!r}"
            ) from exc
        symbols.append(token)

    if not symbols:
        raise ValueError(f"{path}: no atoms found after charge-spin line")

    return symbols, coordinates


def read_xyz(path: Path) -> Structure:
    """
    Read atomic symbols and Cartesian coordinates from an XYZ file.

    Lines from the 3rd onward are parsed as ``token x y z``. If ``token`` is
    an atomic number, it is converted to an element symbol via
    ``rev_elements_num``.

    Parameters
    ----------
    path : Path
        Path to a ``.xyz`` file.

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.

    Raises
    ------
    ValueError
        If no atoms are found or an unknown atomic number appears.
    """
    symbols: list[str] = []
    coordinates: list[CoordTuple] = []

    with path.open(encoding="utf-8") as fh:
        lines = fh.readlines()

    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        token = parts[0]
        if token.isdigit():
            number = int(token)
            try:
                token = rev_elements_num[number]
            except KeyError as exc:
                raise ValueError(
                    f"{path}: unknown atomic number {number}"
                ) from exc
        symbols.append(token)
        coordinates.append((float(parts[1]), float(parts[2]), float(parts[3])))

    if not symbols:
        raise ValueError(f"{path}: no atoms found in XYZ file")

    return symbols, coordinates


def read_mol(path: Path) -> Structure:
    """
    Read atomic symbols and Cartesian coordinates from an MDL mol (V2000) file.

    Parameters
    ----------
    path : Path
        Path to a ``.mol`` file.

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.

    Raises
    ------
    ValueError
        If the file is V3000 or the atom block cannot be parsed.
    """
    with path.open(encoding="utf-8") as fh:
        lines = fh.readlines()

    if len(lines) < 4:
        raise ValueError(f"{path}: mol file too short")

    for line in lines:
        if "V3000" in line.upper():
            raise ValueError(f"{path}: MDL V3000 is not supported (V2000 only)")

    counts = lines[3]
    try:
        n_atoms = int(counts[:3])
    except ValueError as exc:
        raise ValueError(f"{path}: cannot parse atom count from counts line") from exc

    symbols: list[str] = []
    coordinates: list[CoordTuple] = []
    atom_start = 4
    for line in lines[atom_start : atom_start + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: malformed atom line: {line!r}")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        symbol = parts[3]
        symbols.append(symbol)
        coordinates.append((x, y, z))

    if len(symbols) != n_atoms:
        raise ValueError(
            f"{path}: expected {n_atoms} atoms, got {len(symbols)}"
        )

    return symbols, coordinates


def read_mol2(path: Path) -> Structure:
    """
    Read atomic symbols and Cartesian coordinates from a Tripos mol2 file.

    Element symbols are taken from the ``atom_type`` field before the first
    ``.`` (e.g. ``C.ar`` → ``C``).

    Parameters
    ----------
    path : Path
        Path to a ``.mol2`` file.

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.

    Raises
    ------
    ValueError
        If the ``@<TRIPOS>ATOM`` section is missing or empty.
    """
    symbols: list[str] = []
    coordinates: list[CoordTuple] = []

    with path.open(encoding="utf-8") as fh:
        in_atom = False
        for line in fh:
            stripped = line.strip()
            if stripped.upper().startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if in_atom:
                if stripped.upper().startswith("@<TRIPOS>"):
                    break
                if not stripped:
                    continue
                parts = stripped.split()
                if len(parts) < 6:
                    raise ValueError(f"{path}: malformed ATOM line: {stripped!r}")
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                atom_type = parts[5]
                symbol = atom_type.split(".", 1)[0]
                symbols.append(symbol)
                coordinates.append((x, y, z))

    if not symbols:
        raise ValueError(f"{path}: no atoms found in @<TRIPOS>ATOM section")

    return symbols, coordinates


def read_cif(path: Path) -> Structure:
    """
    Extract the first unique molecule from a CIF file as Cartesian coordinates.

    Uses ``CifReader`` and molecule index ``0``. When ``z_value > 1``, a message
    is printed that the first molecule is used.

    Parameters
    ----------
    path : Path
        Path to a ``.cif`` file.

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.
    """
    cif_reader = CifReader(cif_path=str(path))
    if cif_reader.z_value > 1:
        print(
            f"{path}: z_value={cif_reader.z_value}; "
            "using the first unique molecule (index 0)"
        )

    raw_symbols = cif_reader.unique_symbols[0]
    symbols = [str(s) for s in np.asarray(raw_symbols).tolist()]
    cart = cif_reader.convert_frac_to_cart(cif_reader.unique_coords[0])
    coordinates = [
        (float(row[0]), float(row[1]), float(row[2])) for row in np.asarray(cart)
    ]
    return symbols, coordinates


def read_structure(path: Path) -> Structure:
    """
    Dispatch structure reading by file suffix (case-insensitive).

    Parameters
    ----------
    path : Path
        Path to a structure file (``.gjf``, ``.com``, ``.xyz``, ``.mol``,
        ``.mol2``, or ``.cif``).

    Returns
    -------
    tuple[list[str], list[tuple[float, float, float]]]
        ``(symbols, coordinates)`` in Angstrom.

    Raises
    ------
    ValueError
        If the suffix is not supported.
    """
    suffix = path.suffix.lower()
    if suffix in {".gjf", ".com"}:
        return read_gjf(path)
    if suffix == ".xyz":
        return read_xyz(path)
    if suffix == ".mol":
        return read_mol(path)
    if suffix == ".mol2":
        return read_mol2(path)
    if suffix == ".cif":
        return read_cif(path)
    raise ValueError(
        f"unsupported structure suffix {path.suffix!r}; "
        "expected .gjf, .com, .xyz, .mol, .mol2, or .cif"
    )
