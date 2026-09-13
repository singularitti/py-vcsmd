"""File-writing orchestration. Numerical modules never import this module."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import platform
import re
import shutil
import subprocess
from collections.abc import Callable
from contextlib import ExitStack
from dataclasses import dataclass, fields
from datetime import datetime
from enum import Enum
from pathlib import Path

from . import dynamics
from .config import SimulationConfig, prepare
from .geometry import cell_parameters
from .io import save_config
from .io.checkpoint import load_checkpoint, save_checkpoint
from .models import InitialConditions, NumericalModel, Observables, StepResult
from .models import SimulationState as NumericalState
from .simulation import SimulationState


class RunStatus(Enum):
    COMPLETE = "complete"
    FAILED = "failed"


@dataclass(frozen=True)
class RunReport:
    directory: Path
    status: RunStatus
    completed_steps: int
    requested_steps: int
    final_state: SimulationState | None
    error: str | None = None


def _json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def _safe_error(exc: Exception) -> str:
    if isinstance(exc, OSError):
        return f"{type(exc).__name__}: file operation failed"
    return f"{type(exc).__name__}: {str(exc).replace(str(Path.home()), '~')}"


def _new_directory(root: Path, purpose: str, mode: str) -> Path:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    label = re.sub(r"[^a-z0-9]+", "-", purpose.lower()).strip("-") or "simulation"
    stamp = datetime.now().astimezone().strftime("%Y-%m-%d-%H")
    for attempt in range(1, 10000):
        path = root / f"vcsmd-{label}-{mode}-results-{stamp}-attempt-{attempt:02d}"
        try:
            path.mkdir()
        except FileExistsError:
            continue
        (path / "inputs").mkdir()
        (path / "outputs").mkdir()
        (path / "outputs" / "checkpoints").mkdir()
        return path
    raise FileExistsError("Too many run attempts for this purpose and hour")


def _source_snapshot(directory: Path) -> dict[str, object]:
    """Copy implementation inputs and identify them without local path disclosure."""
    package = Path(__file__).parent
    destination = directory / "inputs" / "source" / "vcsmd"
    hashes = {}
    for source in sorted(package.rglob("*.py")):
        relative = source.relative_to(package)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[relative.as_posix()] = hashlib.sha256(source.read_bytes()).hexdigest()
    project = package.parent.parent
    for name in ("pyproject.toml", "uv.lock"):
        candidate = project / name
        if candidate.is_file():
            target = directory / "inputs" / "source" / name
            shutil.copyfile(candidate, target)
            hashes[name] = hashlib.sha256(candidate.read_bytes()).hexdigest()
    revision = None
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=project,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            revision = result.stdout.strip()
    except OSError:
        pass
    versions = {}
    for name in ("vcsmd", "numpy", "pint", "PyYAML", "tomli", "tomli-w"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return {
        "source_provenance": "copied source files of the executing package",
        "source_sha256": hashes,
        "git_revision": revision,
        "python_version": platform.python_version(),
        "dependencies": versions,
    }


def _unit_for(name: str) -> str:
    if name == "step_index":
        return ""
    if name == "time":
        return "rydberg_time"
    if "temperature" in name:
        return "K"
    if name == "volume":
        return "bohr^3"
    if name in ("pressure", "average_pressure"):
        return "Ry/bohr^3"
    return "Ry"


class _SeriesWriter:
    def __init__(self, directory: Path, model: NumericalModel) -> None:
        self.model = model
        self.stack = ExitStack()
        self.observation_names = [item.name for item in fields(Observables)]
        self.observations = self._writer(
            directory / "observables.csv",
            [
                f"{name} [{_unit_for(name)}]" if _unit_for(name) else name
                for name in self.observation_names
            ],
        )
        self.cells = self._writer(
            directory / "cell_history.csv",
            [
                "step_index",
                "time [rydberg_time]",
                *[f"cell_{i}{j} [bohr]" for i in range(3) for j in range(3)],
                *[
                    f"cell_velocity_{i}{j} [bohr/rydberg_time]"
                    for i in range(3)
                    for j in range(3)
                ],
                "a [bohr]",
                "b [bohr]",
                "c [bohr]",
                "alpha [degree]",
                "beta [degree]",
                "gamma [degree]",
            ],
        )
        self.trajectory = self._writer(
            directory / "trajectory.csv",
            [
                "step_index",
                "particle_index",
                "species",
                "fractional_x",
                "fractional_y",
                "fractional_z",
                "x [bohr]",
                "y [bohr]",
                "z [bohr]",
                "peculiar_vx [bohr/rydberg_time]",
                "peculiar_vy [bohr/rydberg_time]",
                "peculiar_vz [bohr/rydberg_time]",
                "fractional_vx [1/rydberg_time]",
                "fractional_vy [1/rydberg_time]",
                "fractional_vz [1/rydberg_time]",
            ],
        )
        self.events = self._writer(
            directory / "events.csv",
            ["step_index", "kind", "scale", "particle_components", "cell_components"],
        )

    def _writer(self, path: Path, columns: list[str]) -> csv.writer:
        handle = self.stack.enter_context(path.open("w", newline="", encoding="utf-8"))
        writer = csv.writer(handle)
        writer.writerow(columns)
        return writer

    def write(self, result: StepResult, *, include_series: bool) -> None:
        for event in result.events:
            self.events.writerow(
                [
                    event.step_index,
                    event.kind.value,
                    "" if event.scale is None else event.scale,
                    event.particle_components,
                    event.cell_components,
                ]
            )
        if not include_series:
            return
        state = result.state
        self.observations.writerow(
            [getattr(result.observables, name) for name in self.observation_names]
        )
        lengths, angles = cell_parameters(state.cell)
        self.cells.writerow(
            [
                state.step_index,
                state.time,
                *state.cell.ravel(),
                *state.cell_velocity.ravel(),
                *lengths,
                *angles,
            ]
        )
        positions = state.fractional_positions @ state.cell.T
        velocities = state.fractional_velocities @ state.cell.T
        for i, species in enumerate(self.model.species):
            self.trajectory.writerow(
                [
                    state.step_index,
                    i,
                    species,
                    *state.fractional_positions[i],
                    *positions[i],
                    *velocities[i],
                    *state.fractional_velocities[i],
                ]
            )

    def close(self) -> None:
        self.stack.close()


def _execute(
    model: NumericalModel,
    *,
    initial: InitialConditions | None,
    state: NumericalState | None,
    steps: int,
    seed: int,
    output_root: Path,
    purpose: str,
    output_interval: int,
    checkpoint_interval: int,
    config: SimulationConfig | None,
    source_file: Path | None,
    progress: Callable[[int, int], None] | None,
) -> RunReport:
    for name, value in (
        ("steps", steps),
        ("output_interval", output_interval),
        ("checkpoint_interval", checkpoint_interval),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    directory = _new_directory(output_root, purpose, model.mode.value)
    output = directory / "outputs"
    initial_step = state.step_index if state is not None else 0
    target_step = initial_step + steps
    metadata = {
        "schema_version": 1,
        "purpose": purpose,
        "created_at": datetime.now().astimezone().isoformat(timespec="hours"),
        "mode": model.mode.value,
        "requested_steps": steps,
        "initial_step": initial_step,
        "target_step": target_step,
        "particles": len(model.masses),
        "seed": seed if initial is not None else None,
        "output_interval": output_interval,
        "checkpoint_interval": checkpoint_interval,
        "status": "running",
        "completed_steps": 0,
        "observables_phase": "before temperature rescaling and quenching",
        "trajectory_phase": "after temperature rescaling and quenching",
        "internal_units": {
            "length": "bohr",
            "time": "rydberg_time",
            "energy": "rydberg",
            "mass": "rydberg_mass",
            "cell_inertia": (
                "rydberg_mass / bohr ** 4"
                if model.mode.metric_cell or model.mode.strain_cell
                else "rydberg_mass"
            )
            if model.mode.variable_cell
            else None,
        },
    }
    writer = None
    error = None
    status = RunStatus.COMPLETE
    completed = 0
    try:
        metadata.update(_source_snapshot(directory))
        if config is not None:
            save_config(config, directory / "inputs" / "configuration.json")
        if source_file is not None:
            source_file = Path(source_file)
            name = (
                "starting_checkpoint.npz"
                if state is not None
                else "original_configuration" + source_file.suffix
            )
            shutil.copyfile(source_file, directory / "inputs" / name)
            metadata["provided_input"] = {
                "copied_as": name,
                "sha256": hashlib.sha256(source_file.read_bytes()).hexdigest(),
            }
        _json(output / "run_metadata.json", metadata)
        (directory / "README.md").write_text(
            f"# VCSMD run\n\nCreated: {metadata['created_at']}\n\n"
            f"Purpose: {purpose}. Mode: {model.mode.value}. "
            f"Requested {steps} steps, from step {initial_step} through {target_step}, "
            f"with {len(model.masses)} particles.\n\n"
            "Inputs are copied under `inputs/`: the normalized configuration or starting checkpoint, "
            "the supplied source file when available, and the Python implementation with dependency metadata. "
            "Source hashes and versions are recorded in `outputs/run_metadata.json`.\n\n"
            "Outputs contain CSV observables, cell history, particle trajectories and controller events; "
            "a complete `checkpoint.npz`; periodic checkpoints under `checkpoints/`; and run metadata. "
            "A failed run preserves its last valid state and records the failure without changing parameters.\n\n"
            "Observables describe each completed integration step BEFORE controller events, matching the "
            "historical measurement order. Trajectory velocities and checkpoints are AFTER those events. "
            "Cell matrices use columns as basis vectors; CSV matrix components are flattened by rows. "
            "Alpha, beta, gamma are the angles between (b,c), (c,a), (a,b). "
            "Positions are wrapped; velocities are peculiar velocities unless labeled otherwise.\n\n"
            "No reduced, trial, or synthetic simulation is created by this execution.\n",
            encoding="utf-8",
        )
        if state is None:
            if initial is None:
                raise ValueError("Initial conditions are required for a new run")
            state = dynamics.initialize(model, initial, seed=seed)
        save_checkpoint(output / "checkpoint.npz", model, state)
        writer = _SeriesWriter(output, model)
        for result in dynamics.simulate(model, state, steps=steps):
            state = result.state
            completed += 1
            writer.write(
                result,
                include_series=state.step_index % output_interval == 0
                or completed == steps,
            )
            if state.step_index % checkpoint_interval == 0:
                save_checkpoint(
                    output / "checkpoints" / f"step-{state.step_index:08d}.npz",
                    model,
                    state,
                )
                save_checkpoint(output / "checkpoint.npz", model, state)
            if progress is not None and (completed % 100 == 0 or completed == steps):
                progress(completed, steps)
        save_checkpoint(output / "checkpoint.npz", model, state)
    except Exception as exc:  # noqa: BLE001 -- execution boundary must preserve failed runs
        status = RunStatus.FAILED
        error = _safe_error(exc)
        if state is not None:
            try:
                save_checkpoint(output / "checkpoint.npz", model, state)
            except (ValueError, OSError) as checkpoint_error:
                metadata["checkpoint_error"] = _safe_error(checkpoint_error)
    finally:
        if writer is not None:
            writer.close()
        metadata.update(
            status=status.value,
            completed_steps=completed,
            error=error,
            final_step=state.step_index if state is not None else initial_step,
        )
        _json(output / "run_metadata.json", metadata)
        with (directory / "README.md").open("a", encoding="utf-8") as handle:
            handle.write(
                f"\nCompletion status: {status.value}. Completed {completed}/{steps} "
                f"requested steps; final simulation step {metadata['final_step']}.\n"
            )
            if error is not None:
                handle.write(f"\nFailure: {error}\n")
    return RunReport(
        directory,
        status,
        completed,
        steps,
        SimulationState(state) if state is not None else None,
        error,
    )


def run(
    config: SimulationConfig,
    *,
    output_root: Path = Path("runs"),
    purpose: str = "simulation",
    source_file: Path | None = None,
    checkpoint_interval: int = 100,
    progress: Callable[[int, int], None] | None = None,
) -> RunReport:
    """Execute all configured steps and write a documented native run folder.

    This application boundary owns file access: it snapshots inputs, streams
    CSV series, and writes periodic and final checkpoints. ``progress`` is an
    optional callback receiving completed and requested step counts. The report
    contains completion status and the public final state; execution failures
    preserve available outputs and the last valid state for inspection.
    """
    model, initial = prepare(config)
    return _execute(
        model,
        initial=initial,
        state=None,
        steps=config.steps,
        seed=config.seed,
        output_root=output_root,
        purpose=purpose,
        output_interval=config.output_interval,
        checkpoint_interval=checkpoint_interval,
        config=config,
        source_file=source_file,
        progress=progress,
    )


def resume(
    checkpoint: Path,
    *,
    steps: int,
    output_root: Path = Path("runs"),
    purpose: str = "continuation",
    output_interval: int = 1,
    checkpoint_interval: int = 100,
    progress: Callable[[int, int], None] | None = None,
) -> RunReport:
    """Continue a native checkpoint for ``steps`` additional integration steps.

    The checkpoint supplies the model, velocities, acceleration histories,
    reference cell, and controller accumulators. No thermal reinitialization
    occurs. Outputs go into a new documented run folder. Bitwise continuation
    assumes the same numerical implementation and software environment.
    """
    model, state = load_checkpoint(Path(checkpoint))
    return _execute(
        model,
        initial=None,
        state=state,
        steps=steps,
        seed=119,
        output_root=output_root,
        purpose=purpose,
        output_interval=output_interval,
        checkpoint_interval=checkpoint_interval,
        config=None,
        source_file=Path(checkpoint),
        progress=progress,
    )
