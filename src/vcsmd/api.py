"""Quantity-bearing public results over the pure normalized solver."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, fields

import numpy as np
import pint

from . import dynamics
from .models import (
    InitialConditions,
    NumericalModel,
    SimulationEvent,
)
from .models import (
    Observables as NumericalObservables,
)
from .models import (
    SimulationState as NumericalState,
)
from .units import (
    ENERGY,
    FRACTIONAL_VELOCITY,
    LENGTH,
    PRESSURE,
    TIME,
    VELOCITY,
    Quantity,
)


@dataclass(frozen=True)
class SimulationState:
    """A read-only state with dimensional fields exposed as Pint quantities.

    ``numerical`` is an explicit escape hatch for advanced numerical work and
    native persistence. Its arrays are independently owned and read-only.
    """

    numerical: NumericalState

    @property
    def cell(self) -> pint.Quantity:
        return Quantity(self.numerical.cell, LENGTH)

    @property
    def fractional_positions(self) -> np.ndarray:
        return self.numerical.fractional_positions

    @property
    def fractional_velocities(self) -> pint.Quantity:
        return Quantity(self.numerical.fractional_velocities, FRACTIONAL_VELOCITY)

    @property
    def positions(self) -> pint.Quantity:
        result = self.numerical.fractional_positions @ self.numerical.cell.T
        result.setflags(write=False)
        return Quantity(result, LENGTH)

    @property
    def peculiar_velocities(self) -> pint.Quantity:
        result = self.numerical.fractional_velocities @ self.numerical.cell.T
        result.setflags(write=False)
        return Quantity(result, VELOCITY)

    @property
    def cartesian_velocities(self) -> pint.Quantity:
        """Full instantaneous Cartesian derivative, including cell deformation."""
        result = (
            self.numerical.fractional_velocities @ self.numerical.cell.T
            + self.numerical.fractional_positions @ self.numerical.cell_velocity.T
        )
        result.setflags(write=False)
        return Quantity(result, VELOCITY)

    @property
    def cell_velocity(self) -> pint.Quantity:
        return Quantity(self.numerical.cell_velocity, VELOCITY)

    @property
    def fractional_acceleration(self) -> pint.Quantity:
        return Quantity(self.numerical.fractional_acceleration, "1 / rydberg_time ** 2")

    @property
    def cell_acceleration(self) -> pint.Quantity:
        return Quantity(self.numerical.cell_acceleration, "bohr / rydberg_time ** 2")

    @property
    def reference_cell(self) -> pint.Quantity:
        return Quantity(self.numerical.reference_cell, LENGTH)

    @property
    def step_index(self) -> int:
        return self.numerical.step_index

    @property
    def time(self) -> pint.Quantity:
        return Quantity(self.numerical.time, TIME)


@dataclass(frozen=True)
class Observables:
    """Measurements for one completed step, before rescaling or quenching.

    Dimensional values are Pint quantities. Atomic temperature uses particle
    kinetic energy; extended temperature also includes cell kinetic energy.
    The running extended temperature is the quantity used by the controller.
    """

    step_index: int
    time: pint.Quantity
    atomic_potential_energy: pint.Quantity
    atomic_kinetic_energy: pint.Quantity
    cell_kinetic_energy: pint.Quantity
    external_pressure_energy: pint.Quantity
    total_energy: pint.Quantity
    pressure: pint.Quantity
    volume: pint.Quantity
    atomic_temperature: pint.Quantity
    extended_temperature: pint.Quantity
    running_extended_temperature: pint.Quantity
    average_potential_energy: pint.Quantity
    average_kinetic_energy: pint.Quantity
    average_pressure: pint.Quantity


def with_units(observables: NumericalObservables) -> Observables:
    """Attach physical units to normalized measurements without changing them."""
    values = {}
    for item in fields(observables):
        name = item.name
        value = getattr(observables, name)
        if name == "step_index":
            values[name] = value
            continue
        if name == "time":
            unit = TIME
        elif name == "volume":
            unit = "bohr ** 3"
        elif "temperature" in name:
            unit = "kelvin"
        elif name in ("pressure", "average_pressure"):
            unit = PRESSURE
        else:
            unit = ENERGY
        values[name] = Quantity(value, unit)
    return Observables(**values)


@dataclass(frozen=True)
class StepResult:
    """An independently owned next state, measurements, and controller events.

    ``observables`` precede controller actions within the step. ``state`` is
    after those actions and is the state to pass to the next step or save.
    Events describe actions for the caller to report; no logging occurs here.
    """

    state: SimulationState
    observables: Observables
    events: tuple[SimulationEvent, ...] = ()


def initialize(
    model: NumericalModel, initial_conditions: InitialConditions, *, seed: int
) -> SimulationState:
    """Create the public state from the two objects returned by ``prepare``.

    ``seed`` is explicit and affects only a local random generator. Thermal
    initialization removes global center-of-mass momentum; provided velocities
    must already satisfy that constraint. The returned state includes Beeman
    history, exposes physical quantities, and owns read-only numerical arrays.
    No input arrays are modified and no files are accessed.
    """
    return SimulationState(dynamics.initialize(model, initial_conditions, seed=seed))


def step(model: NumericalModel, state: SimulationState) -> StepResult:
    """Advance one Beeman step without changing ``model`` or ``state``.

    Pass a state returned by this public API or ``vcsmd.io.load_checkpoint``.
    The result includes dimensional observables and any rescaling or quenching
    events. Its next state can be passed back to this function unchanged.
    """
    result = dynamics.step(model, state.numerical)
    return StepResult(
        SimulationState(result.state), with_units(result.observables), result.events
    )


def simulate(
    model: NumericalModel, initial_state: SimulationState, *, steps: int
) -> Iterator[StepResult]:
    """Yield exactly ``steps`` successive results from the supplied state.

    ``steps`` is a nonnegative count of additional steps. Results are generated
    lazily, so the caller controls iteration and need not retain a trajectory
    in memory. Each result owns its numerical state; retaining a result does
    not cause later steps to modify it. This function performs no file I/O.
    """
    for result in dynamics.simulate(model, initial_state.numerical, steps=steps):
        yield StepResult(
            SimulationState(result.state), with_units(result.observables), result.events
        )
