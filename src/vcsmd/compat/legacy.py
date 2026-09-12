"""Parsing and conversion of the 1991 fixed-format VCSMD files.

This adapter is intentionally one-way.  It recognizes names such as ``inp``
and ``eal`` here, but those names never cross into the numerical package.
The old solver used Fortran ``D`` exponents and fixed-width records; parsing
therefore accepts both the original widths and whitespace-separated records
found in hand-edited or copied runs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PTMASS = 0.5 * 1837.36
_A0 = 0.529177
_RY = 13.6058
_PCV = _A0**3 / _RY / 1.602
_MODE_BY_CODE = {
    "md": "fixed-dyn",
    "mm": "fixed-min",
    "cd": "cell-dyn",
    "cm": "cell-min",
    "nd": "metric-dyn",
    "nm": "metric-min",
    "sd": "strain-dyn",
    "sm": "strain-min",
}
_KNOWN_FILES = (
    "inp",
    "io",
    "e",
    "eal",
    "ave",
    "p",
    "avec",
    "tv",
    "car",
    "out",
    "sip",
)
_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][-+]?\d+)?"
_FLOAT_RE = re.compile(_FLOAT)
_TOKEN_RE = re.compile(rf"(?P<flag>[scthSCTH])|(?P<number>{_FLOAT})")


class LegacyFormatError(ValueError):
    """Raised when a historical input cannot be represented safely."""


@dataclass(frozen=True)
class ConversionReport:
    """Machine-readable account of a legacy conversion.

    ``source`` and the paths in ``native_files``/``source_files`` are useful
    to the application while the conversion is running.  Serializing via
    :meth:`as_dict` emits only relative names for the report itself, avoiding
    machine-specific paths in a run folder.
    """

    source: Path
    destination: Path
    available_data: tuple[str, ...] = ()
    missing_fields: tuple[str, ...] = ()
    malformed_records: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    native_files: tuple[Path, ...] = ()
    source_files: tuple[Path, ...] = ()
    provenance: Mapping[str, Any] = field(default_factory=dict)
    converted_steps: tuple[int, ...] = ()
    completion_status: str = "complete"

    @property
    def problems(self) -> tuple[str, ...]:
        """Missing data, malformed records, and known reliability limitations."""

        return self.missing_fields + self.malformed_records + self.limitations

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-safe report with destination-relative file names."""

        def rel(paths: Iterable[Path]) -> list[str]:
            result: list[str] = []
            for path in paths:
                try:
                    result.append(path.relative_to(self.destination).as_posix())
                except ValueError:
                    result.append(path.name)
            return result

        return {
            "source": self.source.name,
            "destination": self.destination.name,
            "available_data": list(self.available_data),
            "missing_fields": list(self.missing_fields),
            "malformed_records": list(self.malformed_records),
            "limitations": list(self.limitations),
            "native_files": rel(self.native_files),
            "source_files": rel(self.source_files),
            "provenance": dict(self.provenance),
            "converted_steps": list(self.converted_steps),
            "completion_status": self.completion_status,
        }


def _fortran_float(value: str) -> float:
    """Parse a Fortran or ordinary decimal exponent."""

    return float(value.strip().replace("D", "E").replace("d", "e"))


def _without_comment(line: str) -> str:
    # Parenthesized labels are comments in AllInput.txt and in many copied
    # input files.  A title is handled separately and is never passed here.
    return line.split("(", 1)[0]


def _numbers(line: str) -> list[float]:
    return [
        _fortran_float(match.group(0))
        for match in _FLOAT_RE.finditer(_without_comment(line))
    ]


def _integers(line: str) -> list[int]:
    values = _numbers(line)
    if any(value != int(value) for value in values):
        raise LegacyFormatError(f"expected integer fields: {line!r}")
    return [int(value) for value in values]


def _components(line: str, count: int = 3) -> tuple[list[str], list[float]]:
    """Read old ``flag,value`` fields, including the hand-written variants.

    The fixed Fortran format places a one-character operation immediately
    before each value.  Some historical examples put the operation between
    two values (``0.5  s  0.75``); associating a flag with the next value
    accepts both layouts without changing ordinary records.
    """

    flags: list[str] = []
    values: list[float] = []
    pending = ""
    for match in _TOKEN_RE.finditer(_without_comment(line)):
        if match.group("flag") is not None:
            pending = match.group("flag").lower()
            continue
        values.append(_fortran_float(match.group("number")))
        flags.append(pending)
        pending = ""
    if len(values) != count:
        raise LegacyFormatError(f"expected {count} component values: {line!r}")
    return flags, values


def _apply_operation(value: float, flag: str, *, allow_h: bool = True) -> float:
    if not flag:
        return value
    sign = -1.0 if value < 0 else 1.0
    magnitude = abs(value)
    if flag == "s":
        return sign * math.sqrt(magnitude)
    if flag == "c":
        return sign * magnitude ** (1.0 / 3.0)
    if flag == "t":
        return value / 3.0
    if flag == "h" and allow_h:
        # The original reader divides the signed magnitude by sqrt(3); it
        # does not take a square root of the component itself.
        return value / math.sqrt(3.0)
    raise LegacyFormatError(f"unknown coordinate operation {flag!r}")


def _next_nonempty(lines: Sequence[str], index: int) -> tuple[int, str]:
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index >= len(lines):
        raise LegacyFormatError("unexpected end of legacy input")
    return index, lines[index]


def _replicate(
    positions: Sequence[Sequence[float]],
    species: Sequence[str],
    masses: Sequence[float],
    nsc: Sequence[int],
) -> tuple[list[list[float]], list[str], list[float]]:
    if len(nsc) != 3 or any(value <= 0 for value in nsc):
        raise LegacyFormatError(f"supercell dimensions must be positive: {nsc!r}")
    output_positions: list[list[float]] = []
    output_species: list[str] = []
    output_masses: list[float] = []
    for i1 in range(nsc[0]):
        for i2 in range(nsc[1]):
            for i3 in range(nsc[2]):
                for position, atom, mass in zip(positions, species, masses):
                    output_positions.append(
                        [
                            (position[0] + i1) / nsc[0],
                            (position[1] + i2) / nsc[1],
                            (position[2] + i3) / nsc[2],
                        ]
                    )
                    output_species.append(atom)
                    output_masses.append(mass)
    return output_positions, output_species, output_masses


def _parse_legacy_continuation_lines(
    lines: Sequence[str],
    *,
    code: str,
    title: str,
    init_flag: str,
    io_state: Mapping[str, Any],
    start: int,
) -> tuple[dict[str, Any], int]:
    """Parse the input-control tail of an old continuation file.

    In a continuation run the Fortran program reads the cell and atomic state
    from unit 7 (the historical ``io`` file), while consuming placeholder
    records from ``inp``.  Consequently the fresh-run parser cannot safely
    interpret those placeholders as a structure.  We preserve the available
    ``io`` records and mark the imported state partial.
    """

    nonempty = [index for index in range(start, len(lines)) if lines[index].strip()]
    if len(nonempty) < 6:
        raise LegacyFormatError("continuation input has no complete control tail")
    # The tail is stable in the original format: cutoff, image range, step
    # controls, and temperature controls.
    tail_indices = nonempty[-4:]
    tail_values = [_numbers(lines[index]) for index in tail_indices]
    if (
        not tail_values[0]
        or len(tail_values[1]) < 3
        or len(tail_values[2]) < 3
        or len(tail_values[3]) < 3
    ):
        raise LegacyFormatError("continuation input has malformed simulation controls")
    cutoff = tail_values[0][0]
    image_range = [int(value) for value in tail_values[1][:3]]
    steps, interval, max_rescales = (int(value) for value in tail_values[2][:3])
    temperature, tolerance, timestep = tail_values[3][:3]

    # After the two consumed lattice placeholders, cmass/pressure is the next
    # regular input record in crstl.  Keep this intentionally strict: treating
    # an arbitrary placeholder as an inertia would make continuation unsafe.
    pre_tail = [lines[index] for index in nonempty if index < tail_indices[0]]
    if len(pre_tail) < 3:
        raise LegacyFormatError("continuation input has no cell-inertia record")
    mass_pressure = _numbers(pre_tail[2])
    if not mass_pressure:
        raise LegacyFormatError("continuation cell-inertia record is malformed")
    cell_mass = mass_pressure[0]
    pressure = mass_pressure[1] if len(mass_pressure) > 1 else 0.0

    io_cell = io_state.get("cell")
    io_atoms = io_state.get("atoms")
    if (
        not isinstance(io_cell, Mapping)
        or not isinstance(io_cell.get("values"), Sequence)
        or not isinstance(io_atoms, Sequence)
        or not io_atoms
    ):
        raise LegacyFormatError(
            "continuation io does not contain a usable partial structure"
        )
    cell_values = io_cell["values"]
    if len(cell_values) != 3 or any(
        not isinstance(row, Sequence) or len(row) != 3 for row in cell_values
    ):
        raise LegacyFormatError("continuation io has an incomplete cell")
    cell = [[float(value) for value in row] for row in cell_values]
    positions = [
        list(
            atom.get("fractional_position_relative_to_first_atom", {}).get("values", ())
        )
        for atom in io_atoms
    ]
    if any(len(position) != 3 for position in positions):
        raise LegacyFormatError("continuation io has incomplete positions")
    species = [str(atom.get("species", "")) for atom in io_atoms]
    masses = [float(atom["mass"]["value"]) for atom in io_atoms]
    if any(not value for value in species):
        raise LegacyFormatError("continuation io has incomplete species")
    config: dict[str, Any] = {
        "schema_version": 1,
        "mode": _MODE_BY_CODE[code],
        "structure": {
            "cell": {"values": cell, "unit": "bohr"},
            "fractional_positions": positions,
            "masses": {"values": masses, "unit": "rydberg_mass"},
            "species": species,
        },
        "potential": {
            "sigma": {"value": 3.4 / _A0, "unit": "bohr"},
            "epsilon": {"value": 0.0104 / _RY, "unit": "rydberg"},
        },
        "timestep": {"value": timestep, "unit": "rydberg_time"},
        "steps": steps,
        "cutoff": {"value": cutoff, "unit": "bohr"},
        "image_range": image_range,
        "external_pressure": {"value": pressure * _PCV, "unit": "rydberg / bohr ** 3"},
        "cell_inertia": {
            "value": cell_mass * _PTMASS,
            "unit": "rydberg_mass"
            if code in {"md", "mm", "cd", "cm"}
            else "rydberg_mass / bohr ** 4",
        },
        "temperature": {"value": temperature, "unit": "kelvin"},
        "temperature_control": {
            "interval": interval,
            "max_rescales": max_rescales,
            "relative_tolerance": tolerance,
        },
        # ``c`` rethermalizes from the imported geometry.  The old velocity
        # and acceleration columns remain review data in the converter and
        # are deliberately not passed as a runnable native initialization.
        "initialization": {"mode": "thermal", "seed": 119},
        "title": title.strip(),
    }
    return config, tail_indices[-1] + 1


def _parse_legacy_lines(
    lines: Sequence[str],
    start: int = 0,
    *,
    continuation_state: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], int]:
    """Parse one input block and return its config and exclusive end index."""

    index, title_line = _next_nonempty(lines, start)
    title = title_line.rstrip("\r\n")
    index, line = _next_nonempty(lines, index + 1)
    calc_tokens = _without_comment(line).split()
    if not calc_tokens:
        raise LegacyFormatError("missing historical calculation code")
    code = calc_tokens[0].lower()
    if code not in _MODE_BY_CODE:
        raise LegacyFormatError(f"unknown historical calculation code {code!r}")
    index, line = _next_nonempty(lines, index + 1)
    continuation = _without_comment(line).split()
    if len(continuation) < 1:
        raise LegacyFormatError("missing continuation flag")
    init_flag = continuation[0].lower()
    if init_flag not in {"s", "c", "f"}:
        raise LegacyFormatError(f"unknown continuation flag {init_flag!r}")

    if init_flag != "s" and continuation_state is None:
        raise LegacyFormatError(
            "continuation input requires import_legacy_input(path) so its sibling io restart can be read"
        )
    if init_flag != "s":
        return _parse_legacy_continuation_lines(
            lines,
            code=code,
            title=title,
            init_flag=init_flag,
            io_state=continuation_state,
            start=index + 1,
        )

    index, line = _next_nonempty(lines, index + 1)
    scale_values = _numbers(line)
    if not scale_values:
        raise LegacyFormatError("missing lattice scale")
    lattice_scale = scale_values[0]
    index, line = _next_nonempty(lines, index + 1)
    nsc = _integers(line)[:3]
    if len(nsc) != 3:
        raise LegacyFormatError("supercell dimensions require three integers")

    basis_columns: list[list[float]] = []
    for _ in range(3):
        index, line = _next_nonempty(lines, index + 1)
        flags, values = _components(line)
        basis_columns.append(
            [
                lattice_scale * nsc[j] * _apply_operation(value, flags[j])
                for j, value in enumerate(values)
            ]
        )

    index, line = _next_nonempty(lines, index + 1)
    mass_pressure = _numbers(line)
    if not mass_pressure:
        raise LegacyFormatError("missing cell mass")
    cell_mass = mass_pressure[0]
    pressure = mass_pressure[1] if len(mass_pressure) > 1 else 0.0

    index, line = _next_nonempty(lines, index + 1)
    ntype_values = _integers(line)
    if not ntype_values or ntype_values[0] <= 0:
        raise LegacyFormatError("invalid number of atom types")
    ntype = ntype_values[0]
    positions: list[list[float]] = []
    species: list[str] = []
    masses: list[float] = []
    for _ in range(ntype):
        index, line = _next_nonempty(lines, index + 1)
        record = _without_comment(line).split()
        if len(record) < 3:
            raise LegacyFormatError(f"malformed atom-type record: {line!r}")
        try:
            count = int(record[0])
            atomic_mass = _fortran_float(record[-1])
        except (ValueError, IndexError) as exc:
            raise LegacyFormatError(f"malformed atom-type record: {line!r}") from exc
        atom_name = record[1]
        if count <= 0:
            raise LegacyFormatError(f"atom count must be positive: {line!r}")
        for _ in range(count):
            index, line = _next_nonempty(lines, index + 1)
            flags, values = _components(line)
            positions.append(
                [
                    _apply_operation(value, flags[j], allow_h=False)
                    for j, value in enumerate(values)
                ]
            )
            species.append(atom_name)
            masses.append(atomic_mass * _PTMASS)

    positions, species, masses = _replicate(positions, species, masses, nsc)

    index, line = _next_nonempty(lines, index + 1)
    cutoff_values = _numbers(line)
    if not cutoff_values:
        raise LegacyFormatError("missing pair-potential cutoff")
    cutoff = cutoff_values[0]
    index, line = _next_nonempty(lines, index + 1)
    image_range = _integers(line)[:3]
    if len(image_range) != 3:
        raise LegacyFormatError("image range requires three integers")
    index, line = _next_nonempty(lines, index + 1)
    control = _numbers(line)
    if len(control) < 3:
        raise LegacyFormatError("missing simulation control fields")
    steps, interval, max_rescales = (int(control[0]), int(control[1]), int(control[2]))
    index, line = _next_nonempty(lines, index + 1)
    thermal = _numbers(line)
    if len(thermal) < 3:
        raise LegacyFormatError("missing temperature control fields")
    temperature, tolerance, timestep = thermal[:3]

    # Fortran stores lattice vectors as columns.  Input lines are columns;
    # transpose to the native row-major representation of that matrix.
    cell = [[basis_columns[column][row] for column in range(3)] for row in range(3)]
    config: dict[str, Any] = {
        "schema_version": 1,
        "mode": _MODE_BY_CODE[code],
        "structure": {
            "cell": {"values": cell, "unit": "bohr"},
            "fractional_positions": positions,
            "masses": {"values": masses, "unit": "rydberg_mass"},
            "species": species,
        },
        "potential": {
            "sigma": {"value": 3.4 / _A0, "unit": "bohr"},
            "epsilon": {"value": 0.0104 / _RY, "unit": "rydberg"},
        },
        "timestep": {"value": timestep, "unit": "rydberg_time"},
        "steps": steps,
        "cutoff": {"value": cutoff, "unit": "bohr"},
        "image_range": image_range,
        "external_pressure": {"value": pressure * _PCV, "unit": "rydberg / bohr ** 3"},
        "cell_inertia": {
            "value": cell_mass * _PTMASS,
            "unit": "rydberg_mass"
            if code in {"md", "mm", "cd", "cm"}
            else "rydberg_mass / bohr ** 4",
        },
        "temperature": {"value": temperature, "unit": "kelvin"},
        "temperature_control": {
            "interval": interval,
            "max_rescales": max_rescales,
            "relative_tolerance": tolerance,
        },
        "initialization": {"mode": "thermal", "seed": 119},
        "title": title.strip(),
    }
    return config, index + 1


def parse_legacy_input(text: str) -> dict[str, Any]:
    """Parse one historical ``inp`` file into the native configuration dict."""

    if not isinstance(text, str):
        raise TypeError("legacy input must be text")
    lines = text.splitlines()
    config, _ = _parse_legacy_lines(lines)
    return config


def import_legacy_input(path: Path) -> dict[str, Any]:
    """Read a historical input file and return its normalized native mapping.

    A ``c`` continuation imports geometry from its sibling ``io`` file and
    starts with thermal initialization.  A ``f`` follow input is refused: its
    historical state cannot be represented as a complete native checkpoint.
    """

    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    nonempty = [line for line in lines if line.strip()]
    if len(nonempty) < 3:
        raise LegacyFormatError("legacy input is missing required header records")
    flags = _without_comment(nonempty[2]).split()
    init_flag = flags[0].lower() if flags else ""
    if init_flag == "f":
        raise LegacyFormatError(
            "historical follow input cannot become an exact native checkpoint"
        )
    if init_flag != "c":
        return parse_legacy_input(text)
    io_path = path.with_name("io")
    if not io_path.is_file():
        raise LegacyFormatError(
            "historical continuation input requires a sibling io restart file"
        )
    try:
        io_state = _parse_legacy_io(
            io_path.read_text(encoding="utf-8", errors="replace")
        )
    except (OSError, LegacyFormatError) as exc:
        raise LegacyFormatError(
            "historical continuation io restart is malformed"
        ) from exc
    config, _ = _parse_legacy_lines(lines, continuation_state=io_state)
    return config


def split_legacy_examples(text: str) -> dict[int, str]:
    """Extract complete original examples from an ``AllInput.txt`` document.

    The returned values preserve the original lines exactly through the final
    temperature/tolerance/timestep record.  Separators and any later dormant
    tables are excluded.  Blocks that do not parse as complete inputs are
    skipped rather than being presented as runnable configurations.
    """

    lines = text.splitlines(keepends=True)
    starts: list[tuple[int, int]] = []
    marker = re.compile(r"^\s*Input\s+(\d+)\s*:\s*\r?\n?$", re.IGNORECASE)
    for index, line in enumerate(lines):
        match = marker.match(line)
        if match:
            starts.append((index, int(match.group(1))))
    result: dict[int, str] = {}
    for position, (header, number) in enumerate(starts):
        end_limit = (
            starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        )
        block = lines[header + 1 : end_limit]
        try:
            _, end = _parse_legacy_lines([line.rstrip("\r\n") for line in block])
        except LegacyFormatError:
            continue
        # ``end`` includes the final control line and excludes following
        # separators.  Keep a trailing newline when it was present there.
        result[number] = "".join(block[:end])
    return result


def _parse_legacy_io(text: str) -> dict[str, Any]:
    """Parse the final-state portion of historical unit 7 (``io``).

    The routine intentionally labels positions and velocity columns as
    partial: positions are written relative to the first atom and old builds
    contain duplicated lattice derivative columns.  This data is useful for
    review and initialization diagnostics, but it is never an exact restart.
    """

    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 11:
        raise LegacyFormatError("legacy io is too short to contain a final state")
    scale = _numbers(lines[0])
    if not scale:
        raise LegacyFormatError("legacy io has no scale record")
    cursor = 1
    matrices: list[list[float]] = []
    for _ in range(9):
        values = _numbers(lines[cursor])
        if len(values) < 3:
            raise LegacyFormatError("legacy io has malformed cell state")
        matrices.append(values[:3])
        cursor += 1
    try:
        type_fields = _integers(lines[cursor])
        ntype = type_fields[0]
    except (LegacyFormatError, IndexError) as exc:
        raise LegacyFormatError("legacy io has no atom type count") from exc
    if ntype <= 0:
        raise LegacyFormatError("legacy io has an invalid atom type count")
    cursor += 1
    atoms: list[dict[str, Any]] = []
    for _ in range(ntype):
        if cursor >= len(lines):
            raise LegacyFormatError("legacy io ends before all atom types")
        fields = _without_comment(lines[cursor]).split()
        if len(fields) < 3:
            raise LegacyFormatError("legacy io has malformed atom type")
        try:
            count = int(fields[0])
            mass = _fortran_float(fields[-1]) * _PTMASS
        except (ValueError, TypeError) as exc:
            raise LegacyFormatError(
                "legacy io has malformed atom count or mass"
            ) from exc
        if count <= 0:
            raise LegacyFormatError("legacy io has an invalid atom count")
        if not math.isfinite(mass) or mass <= 0:
            raise LegacyFormatError("legacy io has an invalid atom mass")
        name = fields[1]
        cursor += 1
        for _ in range(count):
            if cursor + 2 >= len(lines):
                raise LegacyFormatError("legacy io ends in an atom state")
            position = _numbers(lines[cursor])[:3]
            velocity = _numbers(lines[cursor + 1])[:3]
            acceleration = _numbers(lines[cursor + 2])[:3]
            if len(position) != 3 or len(velocity) != 3 or len(acceleration) != 3:
                raise LegacyFormatError("legacy io has malformed atom state columns")
            atoms.append(
                {
                    "species": name,
                    "mass": {"value": mass, "unit": "rydberg_mass"},
                    "fractional_position_relative_to_first_atom": {
                        "values": position,
                        "unit": "dimensionless",
                    },
                    "fractional_rate": {
                        "values": velocity,
                        "unit": "1 / rydberg_time",
                    },
                    "fractional_acceleration_history": {
                        "values": acceleration,
                        "unit": "1 / rydberg_time ** 2",
                    },
                }
            )
            cursor += 3
    columns = matrices[:3]
    cell = [
        [scale[0] * columns[column][row] for column in range(3)] for row in range(3)
    ]
    cell_velocity = [
        [matrices[3 + column][row] for column in range(3)] for row in range(3)
    ]
    cell_acceleration = [
        [matrices[6 + column][row] for column in range(3)] for row in range(3)
    ]
    return {
        "scale": {"value": scale[0], "unit": "bohr"},
        "cell": {"values": cell, "unit": "bohr"},
        "raw_cell_records": matrices,
        "cell_velocity": {
            "values": cell_velocity,
            "unit": "bohr / rydberg_time",
        },
        "cell_acceleration_history": {
            "values": cell_acceleration,
            "unit": "bohr / rydberg_time ** 2",
        },
        "atoms": atoms,
        "complete": False,
        "state_quality": {
            "complete": False,
            "missing_fields": [
                "absolute fractional origin",
                "reference-cell geometry",
                "previous acceleration required by Beeman integration",
                "temperature-controller counters and averaging accumulators",
            ],
            "damaged_fields": [
                "historical io contains duplicated cell-derivative columns in the source writer"
            ],
        },
    }


def _fixed_output_record(line: str, value_count: int) -> tuple[list[float], int] | None:
    """Parse one of the solver's ``d12.5,...,i6`` records by field width."""

    width = 12 * value_count + 6
    for offset in (1, 0):
        if len(line) < offset + width:
            continue
        fields = [
            line[offset + 12 * index : offset + 12 * (index + 1)]
            for index in range(value_count)
        ]
        step_field = line[offset + 12 * value_count : offset + width]
        if line[offset + width :].strip():
            continue
        try:
            values = [_fortran_float(field) for field in fields]
            step = int(step_field.strip())
        except (ValueError, TypeError):
            continue
        if not all(math.isfinite(value) for value in values):
            continue
        return values, step
    return None


def _car_header_is_valid(lines: Sequence[str]) -> bool:
    """Recognize the three-line ``car`` header in either historical order."""

    if len(lines) < 3:
        return False
    try:
        first, second, third = (_integers(lines[index]) for index in range(3))
    except LegacyFormatError:
        return False

    # The active writer emits nstep, ntype, and the per-type atom counts.
    if (
        len(first) == 1
        and first[0] > 0
        and len(second) == 1
        and second[0] > 0
        and len(third) == second[0]
        and all(value > 0 for value in third)
    ):
        return True
    # Some documentation describes the same header as ntype, counts, nstep.
    if len(first) == 1 and first[0] > 0 and len(second) == first[0] and len(third) == 1:
        return third[0] > 0 and all(value > 0 for value in second)
    return False


def _read_output_records(
    path: Path, value_count: int
) -> tuple[dict[int, list[float]], list[str]]:
    records: dict[int, list[float]] = {}
    malformed: list[str] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return records, [f"{path.name}: unable to read legacy output"]
    header_pending = True
    for line_number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        if header_pending:
            # All output tapes begin with a header containing nstep.
            header_pending = False
            try:
                header_values = _integers(line)
            except LegacyFormatError:
                header_values = []
            if len(header_values) != 1 or header_values[0] <= 0:
                malformed.append(f"{path.name}:{line_number}: malformed output header")
            continue
        fixed = _fixed_output_record(line, value_count)
        if fixed is not None:
            values, step = fixed
            if step <= 0:
                malformed.append(
                    f"{path.name}:{line_number}: invalid nonpositive step {step}"
                )
            elif step in records:
                malformed.append(
                    f"{path.name}:{line_number}: duplicate step {step}; kept first record"
                )
            else:
                records[step] = values
            continue
        # Copied outputs are sometimes whitespace-separated.  Require an
        # explicit integer final token, so a malformed packed record cannot
        # silently turn an exponent or trailing numeric field into a step.
        tokens = _without_comment(line).split()
        if len(tokens) == value_count + 1 and re.fullmatch(r"[+-]?\d+", tokens[-1]):
            try:
                values = [_fortran_float(token) for token in tokens[:-1]]
                step = int(tokens[-1])
            except (ValueError, TypeError):
                values = []
                step = -1
            if len(values) == value_count and all(
                math.isfinite(value) for value in values
            ):
                if step <= 0:
                    malformed.append(
                        f"{path.name}:{line_number}: invalid nonpositive step {step}"
                    )
                elif step in records:
                    malformed.append(
                        f"{path.name}:{line_number}: duplicate step {step}; kept first record"
                    )
                else:
                    records[step] = values
                continue
        malformed.append(f"{path.name}:{line_number}: malformed fixed-width record")
    return records, malformed


def _write_csv(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fieldnames), extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _legacy_initialization_flag(path: Path) -> str:
    """Read only the continuation marker from a historical input file."""

    try:
        lines = [
            line
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip()
        ]
    except OSError:
        return ""
    if len(lines) < 3:
        return ""
    fields = _without_comment(lines[2]).split()
    return fields[0].lower() if fields else ""


def convert_legacy_run(source: Path, destination: Path) -> ConversionReport:
    """Convert available historical outputs into native JSON/CSV datasets.

    Missing tapes remain missing and malformed rows are listed in the report.
    In particular, a header-only ``car`` file produces a missing trajectory
    entry; no positions or cell vectors are invented from other observables.
    """

    source = Path(source)
    destination = Path(destination)
    if not source.is_dir():
        raise NotADirectoryError(source)
    source_resolved = source.resolve()
    destination_resolved = destination.resolve()
    if destination_resolved == source_resolved:
        raise ValueError("conversion destination must differ from source")
    try:
        destination_resolved.relative_to(source_resolved)
    except ValueError:
        pass
    else:
        raise ValueError("conversion destination must not be inside source")
    if destination.exists() and (
        not destination.is_dir() or any(destination.iterdir())
    ):
        raise FileExistsError("conversion destination must be an empty directory")
    destination.mkdir(parents=True, exist_ok=True)
    inputs = destination / "inputs"
    outputs = destination / "outputs"
    inputs.mkdir()
    outputs.mkdir()
    copied = inputs / "legacy_source"
    copied.mkdir(exist_ok=True)
    source_files: list[Path] = []
    provenance: dict[str, Any] = {}
    for filename in _KNOWN_FILES:
        path = source / filename
        if not path.is_file():
            continue
        target = copied / filename
        shutil.copyfile(path, target)
        source_files.append(path)
        provenance[filename] = {"sha256": _sha256(path), "size": path.stat().st_size}

    malformed: list[str] = []
    available: list[str] = []
    missing: list[str] = []
    limitations: list[str] = []
    native_files: list[Path] = []
    if (source / "sip").is_file():
        available.append(
            "sip (archival only; excluded from computational configuration)"
        )
    input_path = source / "inp"
    if input_path.is_file():
        try:
            native_config = import_legacy_input(input_path)
        except (LegacyFormatError, OSError):
            malformed.append("inp: malformed or unreadable legacy input")
            missing.append("native configuration")
        else:
            # Validate at the application boundary before making the native
            # configuration visible.  This import is intentionally local so
            # parsing itself remains independent of the numerical core.
            try:
                from ..io import config_from_mapping

                config_from_mapping(native_config)
            except (TypeError, ValueError):
                malformed.append(
                    "inp: converted configuration failed native schema validation"
                )
                native_config = None
            if native_config is not None:
                configuration_path = outputs / "configuration.json"
                configuration_path.write_text(
                    json.dumps(native_config, indent=2) + "\n", encoding="utf-8"
                )
                native_files.append(configuration_path)
                available.append("configuration")
            if _legacy_initialization_flag(input_path) == "c":
                limitations.append(
                    "continuation geometry was imported with thermal initialization; io is not an exact checkpoint"
                )
    else:
        missing.append("native configuration (inp unavailable)")
    tapes = {
        "e": ("potential_energy", "total_kinetic_energy", "total_energy", "pv", 4),
        "eal": (
            "atomic_potential_energy",
            "atomic_kinetic_energy",
            "atomic_total_energy",
            "cell_potential_energy",
            "cell_kinetic_energy",
            "cell_total_energy",
            6,
        ),
        "ave": ("average_potential_energy", "average_kinetic_energy", 2),
        "p": ("external_pressure", "internal_pressure", "average_pressure", 3),
        "avec": (
            "cell_modulus_1",
            "cell_modulus_2",
            "cell_modulus_3",
            "cell_angle_12",
            "cell_angle_23",
            "cell_angle_31",
            6,
        ),
        "tv": ("cell_volume", "running_extended_temperature", 2),
    }
    data: dict[str, dict[int, list[float]]] = {}
    for filename, description in tapes.items():
        path = source / filename
        if not path.is_file():
            missing.append(f"{filename}: output file unavailable")
            continue
        *names, count = description
        records, errors = _read_output_records(path, int(count))
        data[filename] = records
        malformed.extend(errors)
        if records:
            available.append(filename)
        elif not errors:
            missing.append(f"{filename} contains no data records")

    all_steps = sorted({step for records in data.values() for step in records})
    if all_steps:
        fields = [
            "step_index",
            "potential_energy [rydberg]",
            "total_kinetic_energy [rydberg]",
            "total_energy [rydberg]",
            "pv [rydberg]",
            "atomic_potential_energy [rydberg]",
            "atomic_kinetic_energy [rydberg]",
            "atomic_total_energy [rydberg]",
            "cell_potential_energy [rydberg]",
            "cell_kinetic_energy [rydberg]",
            "cell_total_energy [rydberg]",
            "average_potential_energy [rydberg]",
            "average_kinetic_energy [rydberg]",
            "external_pressure [rydberg / bohr ** 3]",
            "internal_pressure [rydberg / bohr ** 3]",
            "average_pressure [rydberg / bohr ** 3]",
            "cell_volume [bohr ** 3]",
            "running_extended_temperature [kelvin]",
        ]
        header_by_name = {
            field.split(" [", 1)[0]: field for field in fields if " [" in field
        }
        names_by_tape = {
            "e": [
                header_by_name["potential_energy"],
                header_by_name["total_kinetic_energy"],
                header_by_name["total_energy"],
                header_by_name["pv"],
            ],
            "eal": [
                header_by_name["atomic_potential_energy"],
                header_by_name["atomic_kinetic_energy"],
                header_by_name["atomic_total_energy"],
                header_by_name["cell_potential_energy"],
                header_by_name["cell_kinetic_energy"],
                header_by_name["cell_total_energy"],
            ],
            "ave": [
                header_by_name["average_potential_energy"],
                header_by_name["average_kinetic_energy"],
            ],
            "p": [
                header_by_name["external_pressure"],
                header_by_name["internal_pressure"],
                header_by_name["average_pressure"],
            ],
            "tv": [
                header_by_name["cell_volume"],
                header_by_name["running_extended_temperature"],
            ],
        }
        rows: list[dict[str, Any]] = []
        for step in all_steps:
            row: dict[str, Any] = {"step_index": step}
            for filename, names in names_by_tape.items():
                values = data.get(filename, {}).get(step)
                if values is None:
                    continue
                for header, value in zip(names, values):
                    row[header] = value
            rows.append(row)
        observables = outputs / "observables.csv"
        _write_csv(observables, fields, rows)
        native_files.append(observables)
    else:
        missing.append("observables")

    cell_records = data.get("avec", {})
    if cell_records:
        cell_history = outputs / "cell_history.csv"
        _write_csv(
            cell_history,
            [
                "step_index",
                "cell_modulus_1 [bohr]",
                "cell_modulus_2 [bohr]",
                "cell_modulus_3 [bohr]",
                "cell_angle_12 [degree]",
                "cell_angle_23 [degree]",
                "cell_angle_31 [degree]",
            ],
            (
                dict(
                    zip(
                        [
                            "step_index",
                            "cell_modulus_1 [bohr]",
                            "cell_modulus_2 [bohr]",
                            "cell_modulus_3 [bohr]",
                            "cell_angle_12 [degree]",
                            "cell_angle_23 [degree]",
                            "cell_angle_31 [degree]",
                        ],
                        [step, *values],
                    )
                )
                for step, values in sorted(cell_records.items())
            ),
        )
        native_files.append(cell_history)
        limitations.append(
            "avec contains cell moduli and angles only; no cell-vector history was reconstructed"
        )

    car_path = source / "car"
    if car_path.is_file():
        car_lines = [
            line
            for line in car_path.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
            if line.strip()
        ]
        if not _car_header_is_valid(car_lines):
            malformed.append("car: malformed three-line header")
            missing.append("trajectory (car header is malformed)")
            limitations.append(
                "historical car data was not promoted because its header could not be validated"
            )
        elif len(car_lines) == 3:
            missing.append(
                "trajectory (car has only its three-line header; historical positions were not emitted)"
            )
            limitations.append(
                "the historical car writer had its position output commented out"
            )
        else:
            # The old implementation writes a header and atom counts, but no
            # stable per-step position schema. Preserve the raw artifact and
            # explicitly decline to call it a native trajectory.
            missing.append(
                "trajectory (car records are not a complete native trajectory)"
            )
            limitations.append(
                "historical car data was preserved as source bytes, not promoted to trajectory"
            )
    else:
        missing.append("trajectory (car file unavailable)")

    io_path = source / "io"
    if io_path.is_file():
        try:
            restart = _parse_legacy_io(
                io_path.read_text(encoding="utf-8", errors="replace")
            )
        except (LegacyFormatError, OSError):
            malformed.append("io: malformed legacy restart state")
            missing.append("final restart state")
        else:
            restart_path = outputs / "partial_state.json"
            restart_path.write_text(
                json.dumps(restart, indent=2) + "\n", encoding="utf-8"
            )
            native_files.append(restart_path)
            available.append("partial final structure and state (io)")
            limitations.append(
                "io state is partial and cannot serve as an exact native checkpoint"
            )
    else:
        missing.append("final restart state (io unavailable)")

    manifest = inputs / "conversion_manifest.json"
    manifest_payload = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="hours"),
        "purpose": "One-way conversion of historical VCSMD output artifacts",
        "converter_sha256": _sha256(Path(__file__)),
        "source_files": provenance,
    }
    manifest.write_text(json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    native_files.append(manifest)

    status = (
        "complete_with_warnings" if malformed or missing or limitations else "complete"
    )
    report = ConversionReport(
        source=source,
        destination=destination,
        available_data=tuple(dict.fromkeys(available)),
        missing_fields=tuple(dict.fromkeys(missing)),
        malformed_records=tuple(malformed),
        limitations=tuple(dict.fromkeys(limitations)),
        native_files=tuple(native_files),
        source_files=tuple(source_files),
        provenance=provenance,
        converted_steps=tuple(all_steps),
        completion_status=status,
    )
    report_path = outputs / "conversion_report.json"
    report_path.write_text(
        json.dumps(report.as_dict(), indent=2) + "\n", encoding="utf-8"
    )
    native_files.append(report_path)
    readme = destination / "README.md"
    readme.write_text(
        "# Legacy conversion\n\n"
        f"Purpose: one-way import of historical VCSMD files from `{source.name}`.\n"
        f"Created: {manifest_payload['created_at']} (UTC, hour precision).\n\n"
        f"Completion status: {status}. Imported {len(all_steps)} recorded steps.\n\n"
        "Inputs are copied under `inputs/legacy_source/` with SHA-256 provenance in "
        "`inputs/conversion_manifest.json`. Outputs are native CSV/JSON datasets "
        "under `outputs/`. Missing and partial historical data is listed in "
        "`outputs/conversion_report.json`; "
        "a partial restart is never presented as an exact native checkpoint.\n",
        encoding="utf-8",
    )
    native_files.append(readme)
    # The report file was written before the last two files were added; update
    # once so its inventory is complete.
    report = ConversionReport(**{**asdict(report), "native_files": tuple(native_files)})
    report_path.write_text(
        json.dumps(report.as_dict(), indent=2) + "\n", encoding="utf-8"
    )
    return report
