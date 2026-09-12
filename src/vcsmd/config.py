"""Unit-aware configuration and pure normalization, independent of file formats."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pint

from .models import (
    InitialConditions,
    InitializationMode,
    NumericalModel,
    SimulationMode,
    TemperatureControl,
    readonly,
)
from .units import (
    BOLTZMANN,
    CELL_INERTIA,
    ENERGY,
    FRACTIONAL_VELOCITY,
    LENGTH,
    MASS,
    PRESSURE,
    TIME,
    VELOCITY,
    Quantity,
    magnitude,
    scalar,
)


@dataclass(frozen=True)
class Structure:
    """Cell columns and fractional particle positions; masses are per particle."""

    cell: pint.Quantity
    fractional_positions: np.ndarray
    masses: pint.Quantity
    species: tuple[str, ...]

    def __post_init__(self) -> None:
        cell = readonly(magnitude(self.cell, LENGTH, name="cell"), shape=(3, 3))
        positions = readonly(self.fractional_positions)
        if positions.ndim != 2 or positions.shape[1:] != (3,) or not len(positions):
            raise ValueError("fractional_positions must have shape (N, 3)")
        masses = readonly(
            magnitude(self.masses, MASS, name="masses"), shape=(len(positions),)
        )
        species = tuple(self.species)
        if len(species) != len(positions) or any(
            not isinstance(s, str) or not s for s in species
        ):
            raise ValueError("species must provide one nonempty label per particle")
        if np.any(masses <= 0):
            raise ValueError("Particle masses must be positive")
        object.__setattr__(self, "cell", Quantity(cell, LENGTH))
        object.__setattr__(self, "fractional_positions", positions)
        object.__setattr__(self, "masses", Quantity(masses, MASS))
        object.__setattr__(self, "species", species)


@dataclass(frozen=True)
class LennardJones:
    sigma: pint.Quantity = field(default_factory=lambda: Quantity(3.4, "angstrom"))
    epsilon: pint.Quantity = field(
        default_factory=lambda: Quantity(0.0104, "electron_volt")
    )

    @classmethod
    def original_argon(cls) -> LennardJones:
        """The numerical parameter choices of the original active force model."""
        return cls(Quantity(3.4 / 0.529177, LENGTH), Quantity(0.0104 / 13.6058, ENERGY))


@dataclass(frozen=True)
class SimulationConfig:
    structure: Structure
    mode: SimulationMode
    timestep: pint.Quantity
    steps: int
    cutoff: pint.Quantity
    image_range: tuple[int, int, int]
    temperature: pint.Quantity = field(default_factory=lambda: Quantity(0.0, "kelvin"))
    external_pressure: pint.Quantity = field(
        default_factory=lambda: Quantity(0.0, PRESSURE)
    )
    cell_inertia: pint.Quantity | None = None
    potential: LennardJones = field(default_factory=LennardJones.original_argon)
    temperature_control: TemperatureControl = field(default_factory=TemperatureControl)
    initialization: InitializationMode = InitializationMode.THERMAL
    fractional_velocities: pint.Quantity | None = None
    cell_velocity: pint.Quantity | None = None
    seed: int = 119
    output_interval: int = 1
    title: str = "VCSMD simulation"

    def __post_init__(self) -> None:
        if not isinstance(self.structure, Structure):
            raise TypeError("structure must be a Structure")
        if not isinstance(self.mode, SimulationMode):
            raise TypeError("mode must be a SimulationMode enum member")
        if not isinstance(self.initialization, InitializationMode):
            raise TypeError("initialization must be an InitializationMode enum member")
        for name in ("steps", "output_interval"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, int)
            or self.seed < 0
        ):
            raise ValueError("seed must be a nonnegative integer")
        if not isinstance(self.title, str):
            raise TypeError("title must be text")
        images = np.asarray(self.image_range)
        if images.shape != (3,) or images.dtype.kind not in "iu" or np.any(images < 0):
            raise ValueError("image_range must contain three nonnegative integers")
        object.__setattr__(self, "image_range", tuple(int(n) for n in images))


def prepare(config: SimulationConfig) -> tuple[NumericalModel, InitialConditions]:
    """Normalize quantities into an explicit immutable numerical model and state input."""
    inertia_unit = (
        CELL_INERTIA if config.mode.metric_cell or config.mode.strain_cell else MASS
    )
    if config.cell_inertia is None:
        if config.mode.variable_cell:
            raise ValueError(
                f"{config.mode.value} requires cell_inertia in {inertia_unit}"
            )
        inertia = 1.0  # Unused by fixed-cell equations; no dimensional guess is made.
    else:
        inertia = scalar(config.cell_inertia, inertia_unit, name="cell_inertia")
    model = NumericalModel(
        mode=config.mode,
        masses=magnitude(config.structure.masses, MASS, name="masses"),
        species=config.structure.species,
        timestep=scalar(config.timestep, TIME, name="timestep"),
        cutoff=scalar(config.cutoff, LENGTH, name="cutoff"),
        image_range=config.image_range,
        sigma=scalar(config.potential.sigma, LENGTH, name="sigma"),
        epsilon=scalar(config.potential.epsilon, ENERGY, name="epsilon"),
        external_pressure=scalar(
            config.external_pressure, PRESSURE, name="external_pressure"
        ),
        cell_inertia=inertia,
        temperature=scalar(config.temperature, "kelvin", name="temperature"),
        boltzmann_constant=BOLTZMANN,
        temperature_control=config.temperature_control,
    )
    initial = InitialConditions(
        cell=magnitude(config.structure.cell, LENGTH, name="cell"),
        fractional_positions=config.structure.fractional_positions,
        mode=config.initialization,
        fractional_velocities=None
        if config.fractional_velocities is None
        else magnitude(
            config.fractional_velocities,
            FRACTIONAL_VELOCITY,
            name="fractional_velocities",
        ),
        cell_velocity=None
        if config.cell_velocity is None
        else magnitude(config.cell_velocity, VELOCITY, name="cell_velocity"),
    )
    return model, initial
