"""Pure initialization and Beeman evolution in normalized numerical units."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from .cell_dynamics import cell_acceleration
from .geometry import cell_geometry, wrap_fractional
from .models import (
    EventKind,
    InitialConditions,
    InitializationMode,
    NumericalModel,
    RunningAverages,
    SimulationEvent,
    SimulationState,
    StepResult,
)
from .observables import measure
from .potentials import ForceResult, evaluate_lennard_jones


def _forces(
    model: NumericalModel, cell: np.ndarray, positions: np.ndarray
) -> ForceResult:
    return evaluate_lennard_jones(
        cell,
        positions,
        sigma=model.sigma,
        epsilon=model.epsilon,
        cutoff=model.cutoff,
        image_range=model.image_range,
    )


def _accelerations(
    model: NumericalModel,
    cell: np.ndarray,
    cell_velocity: np.ndarray,
    fractional_velocities: np.ndarray,
    reference_cell: np.ndarray,
    forces: ForceResult,
) -> tuple[np.ndarray, np.ndarray]:
    acceleration = np.linalg.solve(cell, forces.forces.T).T / model.masses[:, None]
    if not model.mode.variable_cell:
        return acceleration, np.zeros((3, 3), dtype=np.float64)
    geometry = cell_geometry(cell)
    metric_rate = cell_velocity.T @ cell + cell.T @ cell_velocity
    coupling = geometry.inverse_metric @ metric_rate
    acceleration = acceleration - fractional_velocities @ coupling.T
    lattice_acceleration = cell_acceleration(
        model,
        cell,
        cell_velocity,
        reference_cell,
        fractional_velocities @ cell.T,
        forces.virial,
    )
    return acceleration, lattice_acceleration


def initialize(
    model: NumericalModel, initial_conditions: InitialConditions, *, seed: int
) -> SimulationState:
    """Create a deterministic state without consuming external random state.

    The previous acceleration is initialized to the current acceleration. This
    supplies a second-order self-start instead of inventing zero past forces.
    Provided velocities must respect the removed center-of-mass constraint.
    """
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    cell = initial_conditions.cell
    positions = wrap_fractional(initial_conditions.fractional_positions)
    if len(positions) != len(model.masses):
        raise ValueError("Initial particle count and model masses disagree")
    cell_geometry(cell)
    if initial_conditions.mode is InitializationMode.THERMAL:
        velocities = np.zeros_like(positions)
        if model.temperature > 0:
            if len(positions) == 1:
                raise ValueError(
                    "A single atom has no thermal degrees of freedom after center-of-mass removal"
                )
            generator = np.random.default_rng(seed)
            velocities = generator.normal(size=positions.shape) * np.sqrt(
                model.boltzmann_constant * model.temperature / model.masses[:, None]
            )
            velocities -= np.average(velocities, axis=0, weights=model.masses)
            kinetic = 0.5 * np.sum(model.masses[:, None] * velocities**2)
            velocities *= np.sqrt(
                3
                * (len(positions) - 1)
                * model.boltzmann_constant
                * model.temperature
                / (2 * kinetic)
            )
        fractional_velocities = np.linalg.solve(cell, velocities.T).T
    else:
        fractional_velocities = np.array(
            initial_conditions.fractional_velocities, copy=True
        )
        velocities = fractional_velocities @ cell.T
        center_velocity = np.average(velocities, axis=0, weights=model.masses)
        if np.linalg.norm(center_velocity) > 1e-12 * max(
            1.0, float(np.linalg.norm(velocities))
        ):
            raise ValueError(
                "Provided velocities must have zero total center-of-mass momentum"
            )
    cell_velocity = (
        np.zeros((3, 3))
        if initial_conditions.cell_velocity is None
        else np.array(initial_conditions.cell_velocity, copy=True)
    )
    if not model.mode.variable_cell and np.any(cell_velocity != 0):
        raise ValueError("Fixed-cell modes require zero cell velocity")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        forces = _forces(model, cell, positions)
        acceleration, lattice_acceleration = _accelerations(
            model, cell, cell_velocity, fractional_velocities, cell, forces
        )
    return SimulationState(
        cell=cell,
        fractional_positions=positions,
        fractional_velocities=fractional_velocities,
        cell_velocity=cell_velocity,
        fractional_acceleration=acceleration,
        previous_fractional_acceleration=acceleration,
        cell_acceleration=lattice_acceleration,
        previous_cell_acceleration=lattice_acceleration,
        reference_cell=cell,
        rescalings_remaining=model.temperature_control.max_rescales,
    )


def _advance(model: NumericalModel, state: SimulationState) -> StepResult:
    if len(state.fractional_positions) != len(model.masses):
        raise ValueError("State particle count and model masses disagree")
    dt = model.timestep
    displacement = dt * state.fractional_velocities + dt**2 / 6 * (
        4 * state.fractional_acceleration - state.previous_fractional_acceleration
    )
    positions = wrap_fractional(state.fractional_positions + displacement)
    predicted_velocities = (
        state.fractional_velocities + dt * state.fractional_acceleration
    )
    cell = state.cell
    predicted_cell_velocity = state.cell_velocity
    cell_displacement = np.zeros((3, 3))
    if model.mode.variable_cell:
        cell_displacement = dt * state.cell_velocity + dt**2 / 6 * (
            4 * state.cell_acceleration - state.previous_cell_acceleration
        )
        cell = state.cell + cell_displacement
        predicted_cell_velocity = state.cell_velocity + dt * state.cell_acceleration
    forces = _forces(model, cell, positions)
    acceleration, lattice_acceleration = _accelerations(
        model,
        cell,
        predicted_cell_velocity,
        predicted_velocities,
        state.reference_cell,
        forces,
    )
    velocities = displacement / dt + dt / 6 * (
        2 * acceleration + state.fractional_acceleration
    )
    cell_velocity = np.zeros((3, 3))
    if model.mode.variable_cell:
        cell_velocity = cell_displacement / dt + dt / 6 * (
            2 * lattice_acceleration + state.cell_acceleration
        )
    index = state.step_index + 1
    time = state.time + dt
    observation, averages = measure(
        model,
        cell=cell,
        cell_velocity=cell_velocity,
        fractional_velocities=velocities,
        reference_cell=state.reference_cell,
        forces=forces,
        step_index=index,
        time=time,
        previous_averages=state.averages,
    )
    events = []
    remaining = state.rescalings_remaining
    control = model.temperature_control
    # Retain the original running extended-temperature controller. Zero target
    # temperature disables this controller; initialization already sets zero v.
    if model.temperature > 0 and remaining and index % control.interval == 0:
        running_temperature = observation.running_extended_temperature
        relative_error = abs(running_temperature / model.temperature - 1.0)
        if relative_error > control.relative_tolerance:
            scale = (
                np.sqrt(model.temperature / running_temperature)
                if running_temperature > 1e-13
                else 0.0
            )
            velocities = velocities * scale
            if model.mode.variable_cell and not model.mode.minimizes:
                cell_velocity = cell_velocity * scale
            events.append(
                SimulationEvent(
                    EventKind.TEMPERATURE_RESCALED, index, scale=float(scale)
                )
            )
            averages = RunningAverages()
            remaining -= 1
    if model.mode.minimizes:
        particle_mask = acceleration * state.fractional_acceleration < 0
        cell_mask = lattice_acceleration * state.cell_acceleration < 0
        velocities = np.where(particle_mask, 0.0, velocities)
        if model.mode.variable_cell:
            cell_velocity = np.where(cell_mask, 0.0, cell_velocity)
        if np.any(particle_mask) or np.any(cell_mask):
            events.append(
                SimulationEvent(
                    EventKind.VELOCITIES_QUENCHED,
                    index,
                    particle_components=int(np.count_nonzero(particle_mask)),
                    cell_components=int(np.count_nonzero(cell_mask)),
                )
            )
    next_state = SimulationState(
        cell=cell,
        fractional_positions=positions,
        fractional_velocities=velocities,
        cell_velocity=cell_velocity,
        fractional_acceleration=acceleration,
        previous_fractional_acceleration=state.fractional_acceleration,
        cell_acceleration=lattice_acceleration,
        previous_cell_acceleration=state.cell_acceleration,
        reference_cell=state.reference_cell,
        step_index=index,
        time=time,
        averages=averages,
        rescalings_remaining=remaining,
    )
    return StepResult(next_state, observation, tuple(events))


def step(model: NumericalModel, state: SimulationState) -> StepResult:
    """Advance once, returning new arrays and pre-control observables.

    The one-pass velocity-dependent predictor/corrector matches the source's
    method. Stored accelerations are those evaluated at predicted velocities.
    """
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        return _advance(model, state)


def simulate(
    model: NumericalModel, initial_state: SimulationState, *, steps: int
) -> Iterator[StepResult]:
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 0:
        raise ValueError("steps must be a nonnegative integer")
    state = initial_state
    for _ in range(steps):
        result = step(model, state)
        yield result
        state = result.state
