"""The versioned native configuration schema.

All three supported text formats are converted to the same plain mapping
before validation.  Only this module turns serialized mode strings into enum
members; the numerical core consequently never has to interpret file syntax.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

from ..config import LennardJones, SimulationConfig, Structure, prepare
from ..models import InitializationMode, SimulationMode, TemperatureControl
from ..units import parse_quantity, quantity_record

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib  # type: ignore[no-redef]

import tomli_w
import yaml

SCHEMA_VERSION = 1

_ROOT_FIELDS = {
    "schema_version",
    "mode",
    "structure",
    "potential",
    "timestep",
    "steps",
    "cutoff",
    "image_range",
    "external_pressure",
    "cell_inertia",
    "temperature",
    "temperature_control",
    "initialization",
    "output_interval",
    "title",
}
_STRUCTURE_FIELDS = {"cell", "fractional_positions", "masses", "species"}
_POTENTIAL_FIELDS = {"sigma", "epsilon", "preset"}
_CONTROL_FIELDS = {"interval", "max_rescales", "relative_tolerance"}
_INITIALIZATION_FIELDS = {"mode", "seed", "fractional_velocities", "cell_velocity"}


def _mapping(value: Any, *, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _check_fields(value: Mapping[str, Any], allowed: set[str], *, name: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        fields = ", ".join(repr(item) for item in unknown)
        raise ValueError(f"Unknown field(s) in {name}: {fields}")


def _required(value: Mapping[str, Any], field: str, *, name: str) -> Any:
    if field not in value:
        raise ValueError(f"Missing required field {name}.{field}")
    return value[field]


def _mode(value: Any, *, name: str) -> SimulationMode:
    if isinstance(value, SimulationMode):
        return value
    if not isinstance(value, str):
        raise TypeError(f"{name} must be one of the serialized simulation mode values")
    try:
        return SimulationMode(value)
    except ValueError as exc:
        choices = ", ".join(mode.value for mode in SimulationMode)
        raise ValueError(
            f"Invalid {name} {value!r}; expected one of: {choices}"
        ) from exc


def _initialization_mode(value: Any, *, name: str) -> InitializationMode:
    if isinstance(value, InitializationMode):
        return value
    if not isinstance(value, str):
        raise TypeError(f"{name} must be 'thermal' or 'provided'")
    try:
        return InitializationMode(value)
    except ValueError as exc:
        raise ValueError(
            f"Invalid {name} {value!r}; expected 'thermal' or 'provided'"
        ) from exc


def _integer(value: Any, *, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def _fractional_positions(value: Any) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "structure.fractional_positions must be a finite numeric array"
        ) from exc
    if result.ndim != 2 or result.shape[1] != 3 or not len(result):
        raise ValueError("structure.fractional_positions must have shape (N, 3)")
    if not np.all(np.isfinite(result)):
        raise ValueError(
            "structure.fractional_positions must contain only finite values"
        )
    return result


def _quantity(value: Any, *, name: str, array: bool) -> Any:
    """Parse a quantity record while enforcing scalar/array record syntax."""

    record = _mapping(value, name=name)
    expected = {"values", "unit"} if array else {"value", "unit"}
    if set(record) != expected:
        shape = "values, unit" if array else "value, unit"
        raise ValueError(f"{name} requires exactly {{{shape}}}")
    return parse_quantity(record, name=name)


def config_from_mapping(mapping: Mapping[str, Any]) -> SimulationConfig:
    """Parse and validate one native configuration mapping.

    ``mapping`` is copied into domain objects; callers may safely reuse or
    mutate their input mapping after this function returns.  Calling
    :func:`prepare` here makes dimensional and formulation-specific errors
    appear while loading, before execution starts.
    """

    root = _mapping(mapping, name="configuration")
    _check_fields(root, _ROOT_FIELDS, name="configuration")
    version = _required(root, "schema_version", name="configuration")
    if (
        isinstance(version, bool)
        or not isinstance(version, (int, np.integer))
        or int(version) != SCHEMA_VERSION
    ):
        raise ValueError(
            f"Unsupported configuration schema_version {version!r}; expected {SCHEMA_VERSION}"
        )

    structure_map = _mapping(
        _required(root, "structure", name="configuration"), name="structure"
    )
    _check_fields(structure_map, _STRUCTURE_FIELDS, name="structure")
    structure = Structure(
        cell=_quantity(
            _required(structure_map, "cell", name="structure"),
            name="structure.cell",
            array=True,
        ),
        fractional_positions=_fractional_positions(
            _required(structure_map, "fractional_positions", name="structure")
        ),
        masses=_quantity(
            _required(structure_map, "masses", name="structure"),
            name="structure.masses",
            array=True,
        ),
        species=tuple(_required(structure_map, "species", name="structure")),
    )

    potential_map = _mapping(
        _required(root, "potential", name="configuration"), name="potential"
    )
    _check_fields(potential_map, _POTENTIAL_FIELDS, name="potential")
    preset = potential_map.get("preset")
    if preset is not None and preset != "argon":
        raise ValueError("potential.preset must be 'argon' when present")
    has_sigma = "sigma" in potential_map
    has_epsilon = "epsilon" in potential_map
    if preset is not None:
        if has_sigma or has_epsilon:
            raise ValueError(
                "potential.preset cannot be combined with explicit sigma or epsilon"
            )
        potential = LennardJones.original_argon()
    else:
        if not has_sigma or not has_epsilon:
            missing = "sigma" if not has_sigma else "epsilon"
            raise ValueError(f"Missing required field potential.{missing}")
        potential = LennardJones(
            sigma=_quantity(
                potential_map["sigma"], name="potential.sigma", array=False
            ),
            epsilon=_quantity(
                potential_map["epsilon"], name="potential.epsilon", array=False
            ),
        )

    control_map = _mapping(
        root.get("temperature_control", {}), name="temperature_control"
    )
    _check_fields(control_map, _CONTROL_FIELDS, name="temperature_control")
    control = TemperatureControl(
        interval=_integer(
            control_map.get("interval", 100),
            name="temperature_control.interval",
            minimum=1,
        ),
        max_rescales=_integer(
            control_map.get("max_rescales", 0), name="temperature_control.max_rescales"
        ),
        relative_tolerance=float(control_map.get("relative_tolerance", 0.01)),
    )

    initialization_map = _mapping(root.get("initialization", {}), name="initialization")
    _check_fields(initialization_map, _INITIALIZATION_FIELDS, name="initialization")
    initialization = _initialization_mode(
        initialization_map.get("mode", "thermal"), name="initialization.mode"
    )
    fractional_velocities = (
        None
        if "fractional_velocities" not in initialization_map
        else _quantity(
            initialization_map["fractional_velocities"],
            name="initialization.fractional_velocities",
            array=True,
        )
    )
    cell_velocity = (
        None
        if "cell_velocity" not in initialization_map
        else _quantity(
            initialization_map["cell_velocity"],
            name="initialization.cell_velocity",
            array=True,
        )
    )

    image_range_raw = _required(root, "image_range", name="configuration")
    try:
        image_array = np.asarray(image_range_raw)
    except (TypeError, ValueError) as exc:
        raise ValueError("image_range must contain three nonnegative integers") from exc
    if (
        image_array.shape != (3,)
        or image_array.dtype.kind not in "iu"
        or np.any(image_array < 0)
    ):
        raise ValueError("image_range must contain three nonnegative integers")

    config = SimulationConfig(
        structure=structure,
        mode=_mode(_required(root, "mode", name="configuration"), name="mode"),
        timestep=_quantity(
            _required(root, "timestep", name="configuration"),
            name="timestep",
            array=False,
        ),
        steps=_integer(
            _required(root, "steps", name="configuration"), name="steps", minimum=1
        ),
        cutoff=_quantity(
            _required(root, "cutoff", name="configuration"), name="cutoff", array=False
        ),
        image_range=tuple(int(item) for item in image_array),
        temperature=_quantity(
            root.get("temperature", {"value": 0.0, "unit": "kelvin"}),
            name="temperature",
            array=False,
        ),
        external_pressure=_quantity(
            root.get(
                "external_pressure", {"value": 0.0, "unit": "rydberg / bohr ** 3"}
            ),
            name="external_pressure",
            array=False,
        ),
        cell_inertia=(
            None
            if "cell_inertia" not in root
            else _quantity(root["cell_inertia"], name="cell_inertia", array=False)
        ),
        potential=potential,
        temperature_control=control,
        initialization=initialization,
        fractional_velocities=fractional_velocities,
        cell_velocity=cell_velocity,
        seed=_integer(initialization_map.get("seed", 119), name="initialization.seed"),
        output_interval=_integer(
            root.get("output_interval", 1), name="output_interval", minimum=1
        ),
        title=root.get("title", "VCSMD simulation"),
    )
    prepare(config)
    return config


def config_to_mapping(config: SimulationConfig) -> dict[str, Any]:
    """Return the canonical, format-independent representation of ``config``."""

    if not isinstance(config, SimulationConfig):
        raise TypeError("config must be a SimulationConfig")
    prepare(config)
    structure = config.structure
    initialization: dict[str, Any] = {
        "mode": config.initialization.value,
        "seed": config.seed,
    }
    if config.fractional_velocities is not None:
        initialization["fractional_velocities"] = quantity_record(
            config.fractional_velocities
        )
    if config.cell_velocity is not None:
        initialization["cell_velocity"] = quantity_record(config.cell_velocity)
    potential: dict[str, Any] = {
        "sigma": quantity_record(config.potential.sigma),
        "epsilon": quantity_record(config.potential.epsilon),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "mode": config.mode.value,
        "structure": {
            "cell": quantity_record(structure.cell),
            "fractional_positions": np.asarray(structure.fractional_positions).tolist(),
            "masses": quantity_record(structure.masses),
            "species": list(structure.species),
        },
        "potential": potential,
        "timestep": quantity_record(config.timestep),
        "steps": config.steps,
        "cutoff": quantity_record(config.cutoff),
        "image_range": list(config.image_range),
        "external_pressure": quantity_record(config.external_pressure),
        **(
            {"cell_inertia": quantity_record(config.cell_inertia)}
            if config.cell_inertia is not None
            else {}
        ),
        "temperature": quantity_record(config.temperature),
        "temperature_control": {
            "interval": config.temperature_control.interval,
            "max_rescales": config.temperature_control.max_rescales,
            "relative_tolerance": config.temperature_control.relative_tolerance,
        },
        "initialization": initialization,
        "output_interval": config.output_interval,
        "title": config.title,
    }


def load_config(path: Path) -> SimulationConfig:
    """Load JSON, YAML, or TOML based on ``path``'s suffix."""

    path = Path(path)
    suffix = path.suffix.lower()
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError("Unable to read configuration") from exc
    try:
        if suffix == ".json":
            value = json.loads(text)
        elif suffix in {".yaml", ".yml"}:
            value = yaml.safe_load(text)
        elif suffix == ".toml":
            value = tomllib.loads(text)
        else:
            raise ValueError(
                "Configuration format must use .json, .yaml, .yml, or .toml"
            )
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError("Unable to parse configuration") from exc
    return config_from_mapping(value)


def save_config(config: SimulationConfig, path: Path) -> None:
    """Write canonical configuration data in the format selected by suffix."""

    path = Path(path)
    suffix = path.suffix.lower()
    mapping = config_to_mapping(config)
    try:
        if suffix == ".json":
            text = json.dumps(mapping, indent=2, sort_keys=False) + "\n"
        elif suffix in {".yaml", ".yml"}:
            text = yaml.safe_dump(mapping, sort_keys=False)
        elif suffix == ".toml":
            text = tomli_w.dumps(mapping)
        else:
            raise ValueError(
                "Configuration format must use .json, .yaml, .yml, or .toml"
            )
        path.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise ValueError("Unable to write configuration") from exc
