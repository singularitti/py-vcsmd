"""Pure scalar reductions of the generalized mechanical state."""

import numpy as np

from .cell_dynamics import cell_kinetic_energy
from .geometry import cell_geometry
from .models import NumericalModel, Observables, RunningAverages
from .potentials import ForceResult


def temperatures(
    model: NumericalModel,
    atomic_kinetic: float,
    total_kinetic: float,
    average_kinetic: float,
) -> tuple[float, float, float]:
    atomic_dof = 3 * (len(model.masses) - 1)
    extended_dof = atomic_dof + (6 if model.mode.variable_cell else 0)
    atomic = (
        2 * atomic_kinetic / (atomic_dof * model.boltzmann_constant)
        if atomic_dof
        else 0.0
    )
    if extended_dof:
        scale = 2 / (extended_dof * model.boltzmann_constant)
        return atomic, scale * total_kinetic, scale * average_kinetic
    return atomic, 0.0, 0.0


def measure(
    model: NumericalModel,
    *,
    cell: np.ndarray,
    cell_velocity: np.ndarray,
    fractional_velocities: np.ndarray,
    reference_cell: np.ndarray,
    forces: ForceResult,
    step_index: int,
    time: float,
    previous_averages: RunningAverages,
) -> tuple[Observables, RunningAverages]:
    volume = cell_geometry(cell).volume
    peculiar_velocities = fractional_velocities @ cell.T
    atomic_kinetic = float(0.5 * np.sum(model.masses[:, None] * peculiar_velocities**2))
    cell_kinetic = cell_kinetic_energy(model, cell, cell_velocity, reference_cell)
    external_energy = (
        model.external_pressure * volume if model.mode.variable_cell else 0.0
    )
    potential = forces.potential_energy + external_energy
    kinetic = atomic_kinetic + cell_kinetic
    pressure = float((2 * atomic_kinetic + np.trace(forces.virial)) / (3 * volume))
    averages = RunningAverages(
        count=previous_averages.count + 1,
        potential_sum=previous_averages.potential_sum + potential,
        kinetic_sum=previous_averages.kinetic_sum + kinetic,
        pressure_sum=previous_averages.pressure_sum + pressure,
    )
    avg_kinetic = averages.kinetic_sum / averages.count
    atomic_temperature, extended_temperature, running_temperature = temperatures(
        model, atomic_kinetic, kinetic, avg_kinetic
    )
    observation = Observables(
        step_index=step_index,
        time=time,
        atomic_potential_energy=forces.potential_energy,
        atomic_kinetic_energy=atomic_kinetic,
        cell_kinetic_energy=cell_kinetic,
        external_pressure_energy=external_energy,
        total_energy=potential + kinetic,
        pressure=pressure,
        volume=volume,
        atomic_temperature=atomic_temperature,
        extended_temperature=extended_temperature,
        running_extended_temperature=running_temperature,
        average_potential_energy=averages.potential_sum / averages.count,
        average_kinetic_energy=avg_kinetic,
        average_pressure=averages.pressure_sum / averages.count,
    )
    return observation, averages
