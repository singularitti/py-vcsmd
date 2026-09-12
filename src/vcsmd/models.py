"""Immutable data shared by the numerical kernels. All numbers are internal units.

The unit-aware user-facing configuration lives in :mod:`vcsmd.config`.
This module deliberately has no file, configuration-parser, or Pint dependency.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from numbers import Real
from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def finite_scalar(value: Any, *, name: str) -> float:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Real)
        or not np.isfinite(value)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    return float(value)


def readonly(value: Any, *, shape: tuple[int, ...] | None = None) -> FloatArray:
    """Own an immutable copy, including when the source array is already read-only."""
    array = np.array(value, dtype=np.float64, copy=True)
    if shape is not None and array.shape != shape:
        raise ValueError(f"Expected array shape {shape}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError("Numerical state must contain only finite values")
    array.setflags(write=False)
    return array


class SimulationMode(Enum):
    FIXED_DYNAMICS = "fixed-dyn"
    FIXED_MINIMIZATION = "fixed-min"
    CELL_DYNAMICS = "cell-dyn"
    CELL_MINIMIZATION = "cell-min"
    METRIC_DYNAMICS = "metric-dyn"
    METRIC_MINIMIZATION = "metric-min"
    STRAIN_DYNAMICS = "strain-dyn"
    STRAIN_MINIMIZATION = "strain-min"

    @property
    def variable_cell(self) -> bool:
        return self not in (self.FIXED_DYNAMICS, self.FIXED_MINIMIZATION)

    @property
    def minimizes(self) -> bool:
        return self in (
            self.FIXED_MINIMIZATION,
            self.CELL_MINIMIZATION,
            self.METRIC_MINIMIZATION,
            self.STRAIN_MINIMIZATION,
        )

    @property
    def metric_cell(self) -> bool:
        return self in (self.METRIC_DYNAMICS, self.METRIC_MINIMIZATION)

    @property
    def strain_cell(self) -> bool:
        return self in (self.STRAIN_DYNAMICS, self.STRAIN_MINIMIZATION)


class InitializationMode(Enum):
    THERMAL = "thermal"
    PROVIDED = "provided"


class EventKind(Enum):
    TEMPERATURE_RESCALED = "temperature-rescaled"
    VELOCITIES_QUENCHED = "velocities-quenched"


@dataclass(frozen=True)
class TemperatureControl:
    interval: int = 100
    max_rescales: int = 0
    relative_tolerance: float = 0.01

    def __post_init__(self) -> None:
        if (
            isinstance(self.interval, bool)
            or not isinstance(self.interval, int)
            or self.interval <= 0
        ):
            raise ValueError("Temperature check interval must be a positive integer")
        if (
            isinstance(self.max_rescales, bool)
            or not isinstance(self.max_rescales, int)
            or self.max_rescales < 0
        ):
            raise ValueError(
                "Maximum temperature rescalings must be a nonnegative integer"
            )
        tolerance = finite_scalar(self.relative_tolerance, name="relative_tolerance")
        if tolerance < 0:
            raise ValueError("Temperature tolerance must be finite and nonnegative")
        object.__setattr__(self, "relative_tolerance", tolerance)


@dataclass(frozen=True)
class NumericalModel:
    mode: SimulationMode
    masses: FloatArray
    species: tuple[str, ...]
    timestep: float
    cutoff: float
    image_range: tuple[int, int, int]
    sigma: float
    epsilon: float
    external_pressure: float
    cell_inertia: float
    temperature: float
    boltzmann_constant: float
    temperature_control: TemperatureControl = field(default_factory=TemperatureControl)

    def __post_init__(self) -> None:
        if not isinstance(self.mode, SimulationMode):
            raise TypeError("mode must be a SimulationMode enum member")
        masses = readonly(self.masses)
        if masses.ndim != 1 or not len(masses) or np.any(masses <= 0):
            raise ValueError("masses must be a nonempty vector of positive masses")
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "species", tuple(self.species))
        if len(self.species) != len(masses) or any(
            not isinstance(s, str) or not s for s in self.species
        ):
            raise ValueError("One nonempty species label is required per particle")
        images = tuple(self.image_range)
        if len(images) != 3 or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 0
            for n in images
        ):
            raise ValueError("image_range must contain three nonnegative integers")
        object.__setattr__(self, "image_range", tuple(int(n) for n in images))
        for name in (
            "timestep",
            "cutoff",
            "sigma",
            "epsilon",
            "boltzmann_constant",
            "cell_inertia",
        ):
            value = finite_scalar(getattr(self, name), name=name)
            if value <= 0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)
        object.__setattr__(
            self,
            "external_pressure",
            finite_scalar(self.external_pressure, name="external_pressure"),
        )
        temperature = finite_scalar(self.temperature, name="temperature")
        if temperature < 0:
            raise ValueError("temperature must be finite and nonnegative")
        object.__setattr__(self, "temperature", temperature)
        if not isinstance(self.temperature_control, TemperatureControl):
            raise TypeError("temperature_control must be TemperatureControl")


@dataclass(frozen=True)
class InitialConditions:
    cell: FloatArray
    fractional_positions: FloatArray
    mode: InitializationMode = InitializationMode.THERMAL
    fractional_velocities: FloatArray | None = None
    cell_velocity: FloatArray | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell", readonly(self.cell, shape=(3, 3)))
        positions = readonly(self.fractional_positions)
        if positions.ndim != 2 or positions.shape[1] != 3 or not len(positions):
            raise ValueError("fractional_positions must have shape (N, 3)")
        object.__setattr__(self, "fractional_positions", positions)
        if not isinstance(self.mode, InitializationMode):
            raise TypeError("Initialization mode must be an InitializationMode")
        if self.fractional_velocities is not None:
            object.__setattr__(
                self,
                "fractional_velocities",
                readonly(self.fractional_velocities, shape=positions.shape),
            )
        if self.cell_velocity is not None:
            object.__setattr__(
                self, "cell_velocity", readonly(self.cell_velocity, shape=(3, 3))
            )
        if (
            self.mode is InitializationMode.PROVIDED
            and self.fractional_velocities is None
        ):
            raise ValueError("Provided initialization requires fractional velocities")
        if (
            self.mode is InitializationMode.THERMAL
            and self.fractional_velocities is not None
        ):
            raise ValueError(
                "Thermal initialization cannot also provide particle velocities"
            )


@dataclass(frozen=True)
class RunningAverages:
    count: int = 0
    potential_sum: float = 0.0
    kinetic_sum: float = 0.0
    pressure_sum: float = 0.0

    def __post_init__(self) -> None:
        if (
            isinstance(self.count, bool)
            or not isinstance(self.count, int)
            or self.count < 0
        ):
            raise ValueError("Average count must be a nonnegative integer")
        for name in ("potential_sum", "kinetic_sum", "pressure_sum"):
            value = finite_scalar(getattr(self, name), name=name)
            if (name == "kinetic_sum" and value < 0) or (
                self.count == 0 and value != 0
            ):
                raise ValueError("Invalid running energy or pressure accumulator")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class SimulationState:
    cell: FloatArray
    fractional_positions: FloatArray
    fractional_velocities: FloatArray
    cell_velocity: FloatArray
    fractional_acceleration: FloatArray
    previous_fractional_acceleration: FloatArray
    cell_acceleration: FloatArray
    previous_cell_acceleration: FloatArray
    reference_cell: FloatArray
    step_index: int = 0
    time: float = 0.0
    averages: RunningAverages = field(default_factory=RunningAverages)
    rescalings_remaining: int = 0

    def __post_init__(self) -> None:
        positions = readonly(self.fractional_positions)
        if positions.ndim != 2 or positions.shape[1] != 3 or not len(positions):
            raise ValueError("fractional_positions must have shape (N, 3)")
        object.__setattr__(self, "fractional_positions", positions)
        for name in (
            "fractional_velocities",
            "fractional_acceleration",
            "previous_fractional_acceleration",
        ):
            object.__setattr__(
                self, name, readonly(getattr(self, name), shape=positions.shape)
            )
        for name in (
            "cell",
            "cell_velocity",
            "cell_acceleration",
            "previous_cell_acceleration",
            "reference_cell",
        ):
            object.__setattr__(self, name, readonly(getattr(self, name), shape=(3, 3)))
        for name in ("step_index", "rescalings_remaining"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        time = finite_scalar(self.time, name="time")
        if time < 0:
            raise ValueError("time must be nonnegative")
        object.__setattr__(self, "time", time)
        if (
            not isinstance(self.averages, RunningAverages)
            or self.averages.count > self.step_index
        ):
            raise ValueError("Invalid running averages for simulation state")


@dataclass(frozen=True)
class Observables:
    """Internal-unit measurements at a completed step, before control events.

    Use ``vcsmd.with_units`` to obtain their unit-bearing public representation.
    """

    step_index: int
    time: float
    atomic_potential_energy: float
    atomic_kinetic_energy: float
    cell_kinetic_energy: float
    external_pressure_energy: float
    total_energy: float
    pressure: float
    volume: float
    atomic_temperature: float
    extended_temperature: float
    running_extended_temperature: float
    average_potential_energy: float
    average_kinetic_energy: float
    average_pressure: float


@dataclass(frozen=True)
class SimulationEvent:
    kind: EventKind
    step_index: int
    scale: float | None = None
    particle_components: int = 0
    cell_components: int = 0


@dataclass(frozen=True)
class StepResult:
    state: SimulationState
    observables: Observables
    events: tuple[SimulationEvent, ...] = ()
