"""Compare saved final Fortran ``io`` structures with native Python checkpoints.

This is a read-only postprocessing check for the existing full comparison suite.
It reads the original final ``io`` tapes for inputs 1--6 and the copied native
Python ``checkpoint.npz`` files through :func:`vcsmd.load_checkpoint`.  The
Fortran ``io`` writer stores lattice vectors as columns and atom positions
relative to the first atom.  Positions are therefore origin aligned and
compared with minimum-image fractional differences.  Inputs with one atom have
no translation-invariant internal atomic structure and are explicitly marked
uninformative.  Damaged Fortran lattice velocity/acceleration records are not
read or compared.

The script writes ``outputs/final_structures.csv`` and preserves an exact copy
under ``inputs/compare_fortran_structures.py``.  It performs no integration and
creates no additional simulation workload.
"""

import argparse
import csv
import re
import shutil
from pathlib import Path

import numpy as np

from vcsmd import load_checkpoint

FLOAT = re.compile(r"[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[DEde][+-]?\d+)?")

FIELDS = [
    "input",
    "mode",
    "particles",
    "fortran_io",
    "python_checkpoint",
    "python_step",
    "cell_frobenius_relative_error",
    "cell_length_max_error [bohr]",
    "cell_length_max_relative_error",
    "cell_angle_max_error [degree]",
    "position_comparison",
    "position_max_minimum_image_error",
    "position_rms_minimum_image_error",
    "position_mean_minimum_image_error",
    "notes",
]


def numbers(line: str) -> list[float]:
    """Parse Fortran D-exponent numbers without interpreting labels."""

    return [
        float(token.replace("D", "E").replace("d", "e"))
        for token in FLOAT.findall(line)
    ]


def read_io(path: Path) -> tuple[np.ndarray, np.ndarray, int]:
    """Read final cell columns and relative fractional positions from ``io``."""

    lines = path.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError(f"{path}: truncated io file")
    cursor = 0
    scale = numbers(lines[cursor])
    if len(scale) != 1:
        raise ValueError(f"{path}: invalid scale record")
    cursor += 1

    # celq.f writes each lattice vector as one line, with components i=1..3.
    cell = np.empty((3, 3), dtype=np.float64)
    for column in range(3):
        vector = numbers(lines[cursor])
        cursor += 1
        if len(vector) != 3:
            raise ValueError(f"{path}: invalid lattice vector record")
        cell[:, column] = vector

    # Skip the six derivative records.  iio='n' contains historical repeated
    # and permuted columns; only the geometry above is scientifically used.
    cursor += 6
    if cursor >= len(lines):
        raise ValueError(f"{path}: missing atom-type record")
    type_record = numbers(lines[cursor])
    if not type_record:
        raise ValueError(f"{path}: invalid atom-type count")
    ntype = int(type_record[0])
    cursor += 1
    positions: list[list[float]] = []
    for _ in range(ntype):
        if cursor >= len(lines):
            raise ValueError(f"{path}: missing species record")
        atom_record = lines[cursor]
        atom_count_text = atom_record[:5].strip()
        if not atom_count_text:
            raise ValueError(f"{path}: invalid atom count record")
        atom_count = int(atom_count_text)
        cursor += 1
        for _ in range(atom_count):
            if cursor + 2 >= len(lines):
                raise ValueError(f"{path}: truncated atom record")
            position = numbers(lines[cursor])
            if len(position) != 3:
                raise ValueError(f"{path}: invalid fractional position record")
            positions.append(position)
            # Fractional velocity and acceleration do not enter this geometry check.
            cursor += 3
    cell = cell * scale[0]
    positions_array = np.asarray(positions, dtype=np.float64)
    if (
        not np.isfinite(cell).all()
        or not np.isfinite(positions_array).all()
        or np.linalg.det(cell) <= 0
    ):
        raise ValueError("Fortran io contains invalid geometry")
    return cell, positions_array, len(positions)


def cell_lengths(cell: np.ndarray) -> np.ndarray:
    return np.linalg.norm(cell, axis=0)


def cell_angles(cell: np.ndarray) -> np.ndarray:
    """Return Fortran's gamma, alpha, beta ordering in degrees."""

    vectors = [cell[:, i] for i in range(3)]

    def angle(first: np.ndarray, second: np.ndarray) -> float:
        cosine = np.dot(first, second) / np.linalg.norm(first) / np.linalg.norm(second)
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    return np.asarray(
        [
            angle(vectors[0], vectors[1]),
            angle(vectors[1], vectors[2]),
            angle(vectors[2], vectors[0]),
        ],
        dtype=np.float64,
    )


def minimum_image_difference(fortran: np.ndarray, python: np.ndarray) -> np.ndarray:
    """Compare relative fractional positions modulo periodic translations."""

    python_relative = python - python[0]
    difference = python_relative - fortran
    return difference - np.rint(difference)


def suite_rows(suite: Path) -> list[dict[str, str]]:
    with (suite / "outputs" / "run_summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    selected = {int(row["input"]): row for row in rows if 1 <= int(row["input"]) <= 6}
    if sorted(selected) != list(range(1, 7)):
        raise ValueError("comparison suite must contain complete inputs 1 through 6")
    return [selected[number] for number in range(1, 7)]


def compare(suite: Path) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for row in suite_rows(suite):
        number = int(row["input"])
        if row["status"] != "complete":
            raise ValueError(f"input {number}: Fortran run is incomplete")
        fortran_relative = Path(row["fortran_directory"]) / "io"
        fortran_path = suite / "outputs" / fortran_relative
        checkpoint_relative = (
            Path("inputs")
            / "python_runs"
            / f"input-{number:02d}"
            / "outputs"
            / "checkpoint.npz"
        )
        checkpoint_path = suite / checkpoint_relative
        fortran_cell, fortran_positions, fortran_particles = read_io(fortran_path)
        _, python_state = load_checkpoint(checkpoint_path)
        if python_state.step_index != int(row["steps"]):
            raise ValueError(f"input {number}: Python checkpoint is incomplete")
        cell_value = getattr(python_state.cell, "magnitude", python_state.cell)
        python_cell = np.asarray(cell_value, dtype=np.float64)
        python_positions = np.asarray(
            python_state.fractional_positions, dtype=np.float64
        )
        if fortran_particles != len(python_positions) or fortran_particles != int(
            row["particles"]
        ):
            raise ValueError(f"input {number}: Fortran/Python particle count differs")

        cell_error = float(
            np.linalg.norm(python_cell - fortran_cell) / np.linalg.norm(fortran_cell)
        )
        fortran_lengths = cell_lengths(fortran_cell)
        python_lengths = cell_lengths(python_cell)
        length_errors = np.abs(python_lengths - fortran_lengths)
        length_relative = length_errors / fortran_lengths
        angle_error = float(
            np.max(np.abs(cell_angles(python_cell) - cell_angles(fortran_cell)))
        )

        notes: list[str] = []
        if fortran_particles == 1:
            position_comparison = "uninformative_single_atom"
            position_max = position_rms = position_mean = ""
            notes.append("single atom has no translation-invariant internal structure")
        else:
            position_comparison = "origin_aligned_minimum_image"
            differences = minimum_image_difference(fortran_positions, python_positions)
            magnitudes = np.linalg.norm(differences, axis=1)
            position_max = float(np.max(magnitudes))
            position_rms = float(np.sqrt(np.mean(magnitudes**2)))
            position_mean = float(np.mean(magnitudes))
            notes.append("relative positions compared modulo unit-cell translations")
        results.append(
            {
                "input": number,
                "mode": row["mode"],
                "particles": fortran_particles,
                "fortran_io": fortran_relative.as_posix(),
                "python_checkpoint": checkpoint_relative.as_posix(),
                "python_step": int(python_state.step_index),
                "cell_frobenius_relative_error": cell_error,
                "cell_length_max_error [bohr]": float(np.max(length_errors)),
                "cell_length_max_relative_error": float(np.max(length_relative)),
                "cell_angle_max_error [degree]": angle_error,
                "position_comparison": position_comparison,
                "position_max_minimum_image_error": position_max,
                "position_rms_minimum_image_error": position_rms,
                "position_mean_minimum_image_error": position_mean,
                "notes": "; ".join(notes),
            }
        )
    return results


def write_results(suite: Path, rows: list[dict[str, object]]) -> None:
    output = suite / "outputs" / "final_structures.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    destination = suite / "inputs" / "compare_fortran_structures.py"
    if Path(__file__).resolve() != destination.resolve():
        shutil.copyfile(Path(__file__), destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-suite", type=Path, required=True)
    args = parser.parse_args()
    suite = args.comparison_suite
    rows = compare(suite)
    write_results(suite, rows)
    for row in rows:
        print(
            f"input-{int(row['input']):02d} "
            f"cell_rel={float(row['cell_frobenius_relative_error']):.6g} "
            f"length_abs_bohr={float(row['cell_length_max_error [bohr]']):.6g} "
            f"angle_abs_degree={float(row['cell_angle_max_error [degree]']):.6g} "
            f"position={row['position_comparison']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
