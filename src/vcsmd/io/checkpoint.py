"""Exact, non-pickled native restart checkpoints.

The archive deliberately has a small, explicit schema.  Every array is a
named numeric ``.npy`` member of an ``.npz`` archive and the only textual
member is a scalar UTF-8 JSON metadata value.  This prevents an old or
partially compatible archive from being mistaken for a native restart.
"""

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from ..geometry import cell_geometry
from ..models import (
    NumericalModel,
    RunningAverages,
    SimulationMode,
    SimulationState,
    TemperatureControl,
)

_FORMAT = "vcsmd-native-checkpoint"
_VERSION = 1

_MODEL_ARRAYS = {
    "model_masses",
    "model_timestep",
    "model_cutoff",
    "model_sigma",
    "model_epsilon",
    "model_external_pressure",
    "model_cell_inertia",
    "model_temperature",
    "model_boltzmann_constant",
}
_STATE_ARRAYS = {
    "state_cell",
    "state_fractional_positions",
    "state_fractional_velocities",
    "state_cell_velocity",
    "state_fractional_acceleration",
    "state_previous_fractional_acceleration",
    "state_cell_acceleration",
    "state_previous_cell_acceleration",
    "state_reference_cell",
    "state_step_index",
    "state_time",
    "state_average_count",
    "state_average_potential_sum",
    "state_average_kinetic_sum",
    "state_average_pressure_sum",
    "state_rescalings_remaining",
}
_ARRAYS = _MODEL_ARRAYS | _STATE_ARRAYS | {"metadata"}


def _scalar_array(value: Any, *, integer: bool = False) -> np.ndarray:
    dtype = np.int64 if integer else np.float64
    return np.asarray(value, dtype=dtype)


def _metadata(model: NumericalModel) -> str:
    payload = {
        "format": _FORMAT,
        "version": _VERSION,
        "mode": model.mode.value,
        "species": list(model.species),
        "image_range": list(model.image_range),
        "temperature_control": {
            "interval": model.temperature_control.interval,
            "max_rescales": model.temperature_control.max_rescales,
            "relative_tolerance": model.temperature_control.relative_tolerance,
        },
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _archive_arrays(
    model: NumericalModel, state: SimulationState
) -> dict[str, np.ndarray]:
    return {
        "metadata": np.asarray(_metadata(model)),
        "model_masses": np.array(model.masses, dtype=np.float64, copy=True),
        "model_timestep": _scalar_array(model.timestep),
        "model_cutoff": _scalar_array(model.cutoff),
        "model_sigma": _scalar_array(model.sigma),
        "model_epsilon": _scalar_array(model.epsilon),
        "model_external_pressure": _scalar_array(model.external_pressure),
        "model_cell_inertia": _scalar_array(model.cell_inertia),
        "model_temperature": _scalar_array(model.temperature),
        "model_boltzmann_constant": _scalar_array(model.boltzmann_constant),
        "state_cell": np.array(state.cell, dtype=np.float64, copy=True),
        "state_fractional_positions": np.array(
            state.fractional_positions, dtype=np.float64, copy=True
        ),
        "state_fractional_velocities": np.array(
            state.fractional_velocities, dtype=np.float64, copy=True
        ),
        "state_cell_velocity": np.array(
            state.cell_velocity, dtype=np.float64, copy=True
        ),
        "state_fractional_acceleration": np.array(
            state.fractional_acceleration, dtype=np.float64, copy=True
        ),
        "state_previous_fractional_acceleration": np.array(
            state.previous_fractional_acceleration, dtype=np.float64, copy=True
        ),
        "state_cell_acceleration": np.array(
            state.cell_acceleration, dtype=np.float64, copy=True
        ),
        "state_previous_cell_acceleration": np.array(
            state.previous_cell_acceleration, dtype=np.float64, copy=True
        ),
        "state_reference_cell": np.array(
            state.reference_cell, dtype=np.float64, copy=True
        ),
        "state_step_index": _scalar_array(state.step_index, integer=True),
        "state_time": _scalar_array(state.time),
        "state_average_count": _scalar_array(state.averages.count, integer=True),
        "state_average_potential_sum": _scalar_array(state.averages.potential_sum),
        "state_average_kinetic_sum": _scalar_array(state.averages.kinetic_sum),
        "state_average_pressure_sum": _scalar_array(state.averages.pressure_sum),
        "state_rescalings_remaining": _scalar_array(
            state.rescalings_remaining, integer=True
        ),
    }


def _validate_checkpoint_consistency(
    model: NumericalModel, state: SimulationState
) -> None:
    """Validate invariants required for exact native continuation.

    The numerical dataclasses validate local shapes and scalar values.  These
    checks cover relationships between the model and a state that otherwise
    could produce a syntactically valid but non-continuable checkpoint.
    """

    if len(model.masses) != len(state.fractional_positions):
        raise ValueError(
            "Native checkpoint model and state particle counts do not agree"
        )
    try:
        cell_geometry(state.cell)
    except (TypeError, ValueError, np.linalg.LinAlgError) as exc:
        raise ValueError(
            "Native checkpoint current cell is not positive and nonsingular"
        ) from exc
    try:
        cell_geometry(state.reference_cell)
    except (TypeError, ValueError, np.linalg.LinAlgError) as exc:
        raise ValueError(
            "Native checkpoint reference cell is not positive and nonsingular"
        ) from exc

    positions = state.fractional_positions
    if not np.all((positions >= 0.0) & (positions < 1.0)):
        raise ValueError(
            "Native checkpoint fractional positions must be wrapped in [0, 1)"
        )

    expected_time = float(state.step_index) * model.timestep
    if not np.isclose(state.time, expected_time, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError(
            "Native checkpoint time is inconsistent with step_index and timestep"
        )
    if state.averages.count > state.step_index:
        raise ValueError("Native checkpoint average count exceeds step_index")
    if state.rescalings_remaining > model.temperature_control.max_rescales:
        raise ValueError(
            "Native checkpoint rescaling count exceeds the configured limit"
        )

    if not model.mode.variable_cell:
        if np.any(state.cell_velocity != 0.0):
            raise ValueError(
                "Native checkpoint fixed-cell state has nonzero cell velocity"
            )
        if np.any(state.cell_acceleration != 0.0) or np.any(
            state.previous_cell_acceleration != 0.0
        ):
            raise ValueError(
                "Native checkpoint fixed-cell state has nonzero cell acceleration"
            )
        if not np.array_equal(state.cell, state.reference_cell):
            raise ValueError("Native checkpoint fixed-cell state changed its cell")


def save_checkpoint(path: Path, model: NumericalModel, state: SimulationState) -> None:
    """Atomically save a complete native restart checkpoint."""

    if not isinstance(model, NumericalModel):
        raise TypeError("model must be a NumericalModel")
    if not isinstance(state, SimulationState):
        raise TypeError("state must be a SimulationState")
    _validate_checkpoint_consistency(model, state)
    path = Path(path)
    arrays = _archive_arrays(model, state)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=".vcsmd-checkpoint-",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as stream:
            temporary_name = stream.name
            np.savez_compressed(stream, **arrays)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
        temporary_name = None
    except OSError as exc:
        raise ValueError("Unable to write native checkpoint") from exc
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass


def _required_array(archive: Any, name: str) -> np.ndarray:
    if name not in archive:
        raise ValueError(f"Native checkpoint is missing field {name}")
    value = archive[name]
    if not isinstance(value, np.ndarray) or value.dtype.kind not in "fiu":
        raise ValueError(f"Native checkpoint field {name} must be numeric")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"Native checkpoint field {name} contains nonfinite values")
    return value


def _scalar(archive: Any, name: str, *, integer: bool = False) -> int | float:
    value = _required_array(archive, name)
    if value.shape != ():
        raise ValueError(f"Native checkpoint field {name} must be scalar")
    if integer:
        if value.dtype.kind not in "iu":
            raise ValueError(f"Native checkpoint field {name} must be an integer")
        result = int(value)
        if result < 0:
            raise ValueError(f"Native checkpoint field {name} must be nonnegative")
        return result
    if value.dtype.kind not in "fiu":
        raise ValueError(f"Native checkpoint field {name} must be numeric")
    return float(value)


def _metadata_value(archive: Any) -> dict[str, Any]:
    if "metadata" not in archive:
        raise ValueError("Native checkpoint is missing metadata")
    value = archive["metadata"]
    if value.shape != () or value.dtype.kind != "U":
        raise ValueError("Native checkpoint metadata must be a Unicode scalar")
    try:
        metadata = json.loads(str(value.item()))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("Native checkpoint metadata is not valid JSON") from exc
    if not isinstance(metadata, dict):
        raise TypeError("Native checkpoint metadata must be an object")
    if (
        metadata.get("format") != _FORMAT
        or isinstance(metadata.get("version"), bool)
        or metadata.get("version") != _VERSION
    ):
        raise ValueError("Unsupported native checkpoint format or version")
    required = {
        "format",
        "version",
        "mode",
        "species",
        "image_range",
        "temperature_control",
    }
    if set(metadata) != required:
        raise ValueError("Native checkpoint metadata has an invalid schema")
    return metadata


def load_checkpoint(path: Path) -> tuple[NumericalModel, SimulationState]:
    """Load and validate an exact native checkpoint."""

    try:
        archive_context = np.load(Path(path), allow_pickle=False)
    except (OSError, ValueError, TypeError) as exc:
        raise ValueError("Unable to read native checkpoint") from exc
    if not hasattr(archive_context, "files") or not hasattr(archive_context, "close"):
        # ``numpy.load`` also accepts .npy files and returns an ndarray.  Such
        # files are not native checkpoints and must never reach the parser.
        raise ValueError("Unable to read native checkpoint")
    try:
        archive = archive_context
        names = set(archive.files)
        if names != _ARRAYS:
            missing = sorted(_ARRAYS - names)
            unknown = sorted(names - _ARRAYS)
            details: list[str] = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if unknown:
                details.append("unknown " + ", ".join(unknown))
            raise ValueError(
                "Native checkpoint has an invalid field set ("
                + "; ".join(details)
                + ")"
            )
        metadata = _metadata_value(archive)
        try:
            mode = SimulationMode(metadata["mode"])
        except (KeyError, ValueError, TypeError) as exc:
            raise ValueError(
                "Native checkpoint has an invalid simulation mode"
            ) from exc
        species_raw = metadata["species"]
        if not isinstance(species_raw, list) or any(
            not isinstance(item, str) or not item for item in species_raw
        ):
            raise ValueError("Native checkpoint species metadata is invalid")
        image_raw = metadata["image_range"]
        if (
            not isinstance(image_raw, list)
            or len(image_raw) != 3
            or any(
                isinstance(item, bool) or not isinstance(item, int) or item < 0
                for item in image_raw
            )
        ):
            raise ValueError("Native checkpoint image_range metadata is invalid")
        control_raw = metadata["temperature_control"]
        if not isinstance(control_raw, dict) or set(control_raw) != {
            "interval",
            "max_rescales",
            "relative_tolerance",
        }:
            raise ValueError(
                "Native checkpoint temperature control metadata is invalid"
            )
        control = TemperatureControl(
            interval=control_raw["interval"],
            max_rescales=control_raw["max_rescales"],
            relative_tolerance=control_raw["relative_tolerance"],
        )

        masses = _required_array(archive, "model_masses")
        if masses.ndim != 1 or not len(masses) or np.any(masses <= 0):
            raise ValueError("Native checkpoint masses are invalid")
        model = NumericalModel(
            mode=mode,
            masses=masses,
            species=tuple(species_raw),
            timestep=_scalar(archive, "model_timestep"),
            cutoff=_scalar(archive, "model_cutoff"),
            image_range=tuple(image_raw),
            sigma=_scalar(archive, "model_sigma"),
            epsilon=_scalar(archive, "model_epsilon"),
            external_pressure=_scalar(archive, "model_external_pressure"),
            cell_inertia=_scalar(archive, "model_cell_inertia"),
            temperature=_scalar(archive, "model_temperature"),
            boltzmann_constant=_scalar(archive, "model_boltzmann_constant"),
            temperature_control=control,
        )

        positions = _required_array(archive, "state_fractional_positions")
        if (
            positions.ndim != 2
            or positions.shape[1] != 3
            or positions.shape[0] != len(masses)
        ):
            raise ValueError(
                "Native checkpoint particle arrays have inconsistent shapes"
            )
        state = SimulationState(
            cell=_required_array(archive, "state_cell"),
            fractional_positions=positions,
            fractional_velocities=_required_array(
                archive, "state_fractional_velocities"
            ),
            cell_velocity=_required_array(archive, "state_cell_velocity"),
            fractional_acceleration=_required_array(
                archive, "state_fractional_acceleration"
            ),
            previous_fractional_acceleration=_required_array(
                archive, "state_previous_fractional_acceleration"
            ),
            cell_acceleration=_required_array(archive, "state_cell_acceleration"),
            previous_cell_acceleration=_required_array(
                archive, "state_previous_cell_acceleration"
            ),
            reference_cell=_required_array(archive, "state_reference_cell"),
            step_index=_scalar(archive, "state_step_index", integer=True),
            time=_scalar(archive, "state_time"),
            averages=RunningAverages(
                count=_scalar(archive, "state_average_count", integer=True),
                potential_sum=_scalar(archive, "state_average_potential_sum"),
                kinetic_sum=_scalar(archive, "state_average_kinetic_sum"),
                pressure_sum=_scalar(archive, "state_average_pressure_sum"),
            ),
            rescalings_remaining=_scalar(
                archive, "state_rescalings_remaining", integer=True
            ),
        )
        _validate_checkpoint_consistency(model, state)
        return model, state
    except (TypeError, ValueError, OverflowError) as exc:
        if isinstance(exc, ValueError) and str(exc).startswith("Native checkpoint"):
            raise
        raise ValueError("Native checkpoint is malformed") from exc
    finally:
        archive_context.close()
