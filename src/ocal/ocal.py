"""
ocal: single-molecule orbital tail fraction (OTF) calculator.

Assumes a neutral closed-shell molecule (charge/spin fixed to ``0 1``).
Input charge/spin information in structure files is ignored.
Gaussian functional/basis is fixed to B3LYP/6-31G(d,p) for comparability
with previous OTF results; only ``Ocal._create_gjf(method=...)`` may override.
"""
from __future__ import annotations

import argparse
import functools
import os
import subprocess
from datetime import datetime
from pathlib import Path
from time import time

from .utils.cube_reader import calc_otf as _calc_otf_from_cube
from .utils.gaus_log_reader import check_normal_termination as _check_normal_termination
from .utils.gjf_maker import GjfMaker
from .utils.structure_reader import read_structure

print = functools.partial(print, flush=True)

# To make OTF comparable to previous results, the functional and basis set are fixed.
_METHOD = "B3LYP/6-31G(d,p)"
_ORBITALS = ("NHOMO", "HOMO", "LUMO", "NLUMO")
_STRUCTURE_SUFFIXES = {".gjf", ".com", ".xyz", ".mol", ".mol2", ".cif"}
_SUPPORTED_SUFFIXES = _STRUCTURE_SUFFIXES | {".fchk", ".cube"}


def main() -> None:
    """
    Run the ocal CLI: compute OTF for one molecule from a structure, fchk, or cube.
    """
    parser, args = parse_args()

    print("----------------------------------------")
    print(" ocal 0.1.0 (2026/08/06) by Matsui Lab. ")
    print("----------------------------------------")
    print()

    suffix = Path(args.file).suffix.lower()
    if suffix not in _SUPPORTED_SUFFIXES:
        parser.error(
            f"unsupported file extension {Path(args.file).suffix!r}; "
            "expected .gjf, .com, .xyz, .mol, .mol2, .cif, .fchk, or .cube"
        )

    if suffix in {".fchk", ".cube"}:
        if not args.skip_gaussian:
            parser.error(
                "-s/--skip-gaussian is required when input is .fchk or .cube"
            )
    elif args.skip_gaussian:
        parser.error(
            "-s/--skip-gaussian is only valid for .fchk or .cube input"
        )

    print(f"Input File Name: {args.file}")
    Ocal.print_timestamp()
    print()
    start_time = time()

    ocal = Ocal(file=args.file, cpu=args.cpu, mem=args.mem)

    if ocal.suffix in _STRUCTURE_SUFFIXES:
        ocal.create_gjf()
        ocal.run_gaussian()
        ocal.create_fchk()
        ocal.create_cube()
        results = ocal.calc_otf()
    elif ocal.suffix == ".fchk":
        ocal.create_cube()
        results = ocal.calc_otf()
    else:
        results = ocal.calc_otf()

    print()
    print("----------------------------------------")
    print(f"{'orbital':<10} {'OTF':>14} {'density_sum':>14}")
    print("----------------------------------------")
    for label, (otf, density_sum) in results.items():
        print(f"{label:<10} {otf:14.10f} {density_sum:14.10f}")
    print("----------------------------------------")
    print()

    Ocal.print_timestamp()
    _print_elapsed(time() - start_time)


def parse_args() -> tuple[argparse.ArgumentParser, argparse.Namespace]:
    """
    Parse command-line arguments.

    Returns
    -------
    tuple[argparse.ArgumentParser, argparse.Namespace]
        The parser (for ``parser.error``) and parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="ocal",
        description=(
            "Calculate orbital tail fraction (OTF) for a single molecule. "
            "Charge/spin are fixed to 0 1 (neutral singlet)."
        ),
    )
    parser.add_argument("file", help="input file (.gjf/.com/.xyz/.mol/.mol2/.cif/.fchk/.cube)", type=str)
    parser.add_argument("-c", "--cpu", help="setting the number of cpu (default is 4)", type=int, default=4)
    parser.add_argument(
        "-m", "--mem",
        help="setting the number of memory (GB) (default is 10 GB)",
        type=int,
        default=10,
    )
    parser.add_argument(
        "-s", "--skip-gaussian",
        action="store_true",
        help="skip Gaussian calculation (required when input is .fchk or .cube)",
    )

    args = parser.parse_args()
    return parser, args


def _print_elapsed(elapsed_time: float) -> None:
    """
    Print elapsed time in a human-readable form (mcal-style, simplified).

    Parameters
    ----------
    elapsed_time : float
        Elapsed seconds.
    """
    elapsed_time_h = int(elapsed_time // 3600)
    elapsed_time_min = int((elapsed_time - elapsed_time_h * 3600) // 60)
    elapsed_time_sec = int(
        elapsed_time - elapsed_time_h * 3600 - elapsed_time_min * 60
    )
    elapsed_time_ms = (
        elapsed_time
        - elapsed_time_h * 3600
        - elapsed_time_min * 60
        - elapsed_time_sec
    ) * 1000
    if elapsed_time < 1:
        print(f"Elapsed Time: {elapsed_time_ms:.0f} ms")
    elif elapsed_time < 60:
        print(f"Elapsed Time: {elapsed_time_sec} sec")
    elif elapsed_time < 3600:
        print(f"Elapsed Time: {elapsed_time_min} min {elapsed_time_sec} sec")
    else:
        print(
            f"Elapsed Time: {elapsed_time_h} h "
            f"{elapsed_time_min} min {elapsed_time_sec} sec"
        )


def _gaussian_env(base_path: Path) -> dict[str, str]:
    """
    Build environment dict with GAUSS_SCRDIR and GAUSS_ARCDIR.

    Parameters
    ----------
    base_path : Path
        Base path without extension. The parent directory is used
        for both ``GAUSS_SCRDIR`` and ``GAUSS_ARCDIR``.

    Returns
    -------
    dict[str, str]
        Copy of the current environment with Gaussian directories set.
    """
    work_dir = str(base_path.parent)
    env = os.environ.copy()
    env["GAUSS_SCRDIR"] = work_dir
    env["GAUSS_ARCDIR"] = work_dir
    return env


_LOG_SUFFIXES = (".log", ".out")


def _find_gaussian_log(base_path: Path) -> Path | None:
    """
    Locate the Gaussian output file for a job.

    Gaussian on Linux writes ``{base}.log`` while the Windows build writes
    ``{base}.out``; both are searched.

    Parameters
    ----------
    base_path : Path
        Job base path without extension.

    Returns
    -------
    Path | None
        The first existing log path among ``.log`` and ``.out``, or ``None``
        if neither exists.
    """
    for suffix in _LOG_SUFFIXES:
        candidate = Path(f"{base_path}{suffix}")
        if candidate.is_file():
            return candidate
    return None


def _read_alpha_electrons(fchk_path: Path) -> int:
    """
    Read the number of alpha electrons from a formatted checkpoint file.

    Parameters
    ----------
    fchk_path : Path
        Path to a ``.fchk`` file.

    Returns
    -------
    int
        Number of alpha electrons (``na``).

    Raises
    ------
    ValueError
        If the field is not found.
    """
    with fchk_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("Number of alpha electrons"):
                return int(line.split()[-1])
    raise ValueError(f"{fchk_path}: 'Number of alpha electrons' not found")


class Ocal:
    """
    Calculate orbital tail fraction (OTF) for a single molecule.

    Charge and spin multiplicity are fixed to ``0 1`` (neutral singlet).
    Structure-file inputs always run Opt=Tight at B3LYP/6-31G(d,p).
    """

    def __init__(self, file: str, cpu: int = 4, mem: int = 10) -> None:
        """
        Initialize Ocal for one input file.

        Parameters
        ----------
        file : str
            Path to the input file.
        cpu : int
            Number of CPUs for Gaussian / cubegen (default 4).
        mem : int
            Memory in GB for Gaussian (default 10).

        Raises
        ------
        ValueError
            If the file suffix is not supported.
        """
        self.file = Path(file).resolve()
        self.suffix = self.file.suffix.lower()
        self.base_path = self.file.with_suffix("")
        self.cpu = cpu
        self.mem = mem
        self.gau_com = "g16"

        if self.suffix not in _SUPPORTED_SUFFIXES:
            raise ValueError(
                f"unsupported file extension {self.file.suffix!r}; "
                "expected .gjf, .com, .xyz, .mol, .mol2, .cif, .fchk, or .cube"
            )

        # Avoid overwriting input .gjf/.com when export_gjf writes {stem}.gjf
        if self.suffix in {".gjf", ".com"}:
            self.work_base = self.file.parent / f"{self.file.stem}_ocal"
        else:
            self.work_base = self.base_path

    @staticmethod
    def check_normal_termination(log_file: str) -> bool:
        """
        Check whether a Gaussian log ended with normal termination.

        Parameters
        ----------
        log_file : str
            Path to the Gaussian log file.

        Returns
        -------
        bool
            True if the calculation terminated normally.
        """
        return _check_normal_termination(log_file)

    @staticmethod
    def print_timestamp() -> None:
        """Print timestamp."""
        month = {
            1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun",
            7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec",
        }
        dt_now = datetime.now()
        print(
            f"Timestamp: {dt_now.strftime('%a')} {month[dt_now.month]} "
            f"{dt_now.strftime('%d %H:%M:%S %Y')}"
        )

    def create_gjf(self) -> None:
        """
        Read the structure file and write a Gaussian input (``.gjf``).
        """
        symbols, coordinates = read_structure(self.file)
        self._create_gjf(symbols, coordinates)

    def _create_gjf(
        self,
        symbols: list[str],
        coordinates: list[tuple[float, float, float]],
        method: str = _METHOD,
    ) -> None:
        """
        Build and export a Gaussian gjf via GjfMaker.

        Parameters
        ----------
        symbols : list[str]
            Atomic symbols.
        coordinates : list[tuple[float, float, float]]
            Cartesian coordinates in Angstrom.
        method : str
            Functional/basis string (default ``B3LYP/6-31G(d,p)``).
            Kept as an argument for future extension; not exposed on the CLI.
        """
        gjf_maker = GjfMaker()
        gjf_maker.set_resource(self.cpu, self.mem)
        gjf_maker.create_chk_file()
        gjf_maker.opt(tight=True)
        gjf_maker.set_function(method)
        gjf_maker.set_symbols(symbols)
        gjf_maker.set_coordinates(coordinates)
        gjf_maker.set_charge_spin(0, 1)
        gjf_maker.set_title(self.work_base.name)
        gjf_maker.export_gjf(
            file_name=self.work_base.name,
            save_dir=str(self.work_base.parent),
            chk_rwf_name=str(self.work_base),
        )
        print(f"Created: {self.work_base}.gjf")

    def run_gaussian(self, gau_com: str = "g16") -> None:
        """
        Run Gaussian on the generated gjf.

        The subprocess runs with ``cwd`` set to the work directory and a
        relative gjf filename so Windows Gaussian does not concatenate an
        absolute path onto the current directory.

        Parameters
        ----------
        gau_com : str
            Gaussian executable name (default ``g16``). Stored on the instance
            as ``gau_com`` for future CLI extension (e.g. g09).

        Raises
        ------
        RuntimeError
            If the subprocess fails, no ``.log``/``.out`` is found, or the
            log does not show normal termination.
        """
        self.gau_com = gau_com
        work_dir = self.work_base.parent
        gjf_name = f"{self.work_base.name}.gjf"
        command = [gau_com, gjf_name]
        print(f"> {' '.join(command)}")
        try:
            res = subprocess.run(
                command,
                cwd=work_dir,
                env=_gaussian_env(self.work_base),
            )
        except FileNotFoundError as exc:
            raise RuntimeError(f"Failed to execute {' '.join(command)}: {exc}") from exc
        if res.returncode:
            raise RuntimeError(
                f"Failed to execute {' '.join(command)} (returncode={res.returncode})"
            )
        log_path = _find_gaussian_log(self.work_base)
        if log_path is None:
            candidates = ", ".join(f"{self.work_base}{s}" for s in _LOG_SUFFIXES)
            raise RuntimeError(
                f"Gaussian output not found (looked for: {candidates})"
            )
        if not self.check_normal_termination(str(log_path)):
            raise RuntimeError(
                f"Gaussian did not terminate normally: {log_path}"
            )
        print("Gaussian calculation completed")
        print(f" {log_path}")

    def create_fchk(self) -> None:
        """
        Convert the checkpoint file to a formatted checkpoint via formchk.

        Raises
        ------
        RuntimeError
            If formchk fails.
        """
        work_dir = self.work_base.parent
        chk_name = f"{self.work_base.name}.chk"
        fchk_name = f"{self.work_base.name}.fchk"
        command = ["formchk", chk_name, fchk_name]
        print(f"> {' '.join(command)}")
        try:
            res = subprocess.run(
                command,
                cwd=work_dir,
                env=_gaussian_env(self.work_base),
            )
        except FileNotFoundError as exc:
            raise RuntimeError(f"Failed to execute {' '.join(command)}: {exc}") from exc
        if res.returncode:
            raise RuntimeError(
                f"Failed to execute {' '.join(command)} (returncode={res.returncode})"
            )
        print(f"Created: {self.work_base}.fchk")

    def create_cube(self) -> None:
        """
        Generate NHOMO/HOMO/LUMO/NLUMO cube files with cubegen.

        MO indices are derived from the number of alpha electrons ``na`` in
        the fchk: NHOMO=na-1, HOMO=na, LUMO=na+1, NLUMO=na+2.

        Raises
        ------
        RuntimeError
            If cubegen fails for any orbital.
        """
        if self.suffix == ".fchk":
            fchk_src = self.file
        else:
            fchk_src = Path(f"{self.work_base}.fchk")

        na = _read_alpha_electrons(fchk_src)
        mo_indices = {
            "NHOMO": na - 1,
            "HOMO": na,
            "LUMO": na + 1,
            "NLUMO": na + 2,
        }

        work_dir = self.work_base.parent
        fchk_name = fchk_src.name
        for orb in _ORBITALS:
            cube_name = f"{self.work_base.name}_{orb}.cube"
            command = [
                "cubegen",
                str(self.cpu),
                f"MO={mo_indices[orb]}",
                fchk_name,
                cube_name,
                "-2",
                "h",
            ]
            print(f"> {' '.join(command)}")
            try:
                res = subprocess.run(
                    command,
                    cwd=work_dir,
                    env=_gaussian_env(self.work_base),
                )
            except FileNotFoundError as exc:
                raise RuntimeError(
                    f"Failed to execute {' '.join(command)}: {exc}"
                ) from exc
            if res.returncode:
                raise RuntimeError(
                    f"Failed to execute {' '.join(command)} "
                    f"(returncode={res.returncode})"
                )
            print(f"Created: {self.work_base}_{orb}.cube")

    def calc_otf(self) -> dict[str, tuple[float, float]]:
        """
        Compute OTF and density_sum for the target cube file(s).

        Returns
        -------
        dict[str, tuple[float, float]]
            Mapping of orbital label (or cube stem) to ``(otf, density_sum)``.
        """
        results: dict[str, tuple[float, float]] = {}
        if self.suffix == ".cube":
            label = self.file.stem
            results[label] = _calc_otf_from_cube(self.file)
            return results

        for orb in _ORBITALS:
            cube_path = Path(f"{self.work_base}_{orb}.cube")
            results[orb] = _calc_otf_from_cube(cube_path)
        return results


if __name__ == "__main__":
    main()
