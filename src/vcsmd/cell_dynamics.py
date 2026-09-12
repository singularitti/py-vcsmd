"""Cell metrics and generalized forces, expressed as small matrix operations."""

from __future__ import annotations

import numpy as np

from .geometry import cell_geometry
from .models import NumericalModel


def cell_kinetic_energy(
    model: NumericalModel,
    cell: np.ndarray,
    cell_velocity: np.ndarray,
    reference_cell: np.ndarray,
) -> float:
    if not model.mode.variable_cell:
        return 0.0
    weighted = cell_velocity
    if model.mode.metric_cell or model.mode.strain_cell:
        metric_cell = reference_cell if model.mode.strain_cell else cell
        cofactor = cell_geometry(metric_cell).cofactor
        weighted = cell_velocity @ cofactor.T
    return float(0.5 * model.cell_inertia * np.sum(weighted * weighted))


def cell_acceleration(
    model: NumericalModel,
    cell: np.ndarray,
    cell_velocity: np.ndarray,
    reference_cell: np.ndarray,
    peculiar_velocities: np.ndarray,
    virial: np.ndarray,
) -> np.ndarray:
    """Generalized cell acceleration including the common strain projection.

    For metric dynamics F = C.T @ C, where C = det(H) H**(-T).
    The metric derivative term is the gradient of tr(Hdot.T Hdot F)/2
    with Hdot held fixed. Its closed form includes both symmetric terms.
    """
    if not model.mode.variable_cell:
        return np.zeros((3, 3), dtype=np.float64)
    geometry = cell_geometry(cell)
    kinetic_stress = (peculiar_velocities.T * model.masses) @ peculiar_velocities
    stress = (kinetic_stress + virial) / geometry.volume
    driving = (
        (stress - model.external_pressure * np.eye(3))
        @ geometry.cofactor
        / model.cell_inertia
    )
    if model.mode.metric_cell:
        inverse_transpose = geometry.cofactor / geometry.volume
        metric = geometry.cofactor.T @ geometry.cofactor
        inverse_metric = geometry.metric / geometry.volume**2
        velocity_metric = cell_velocity.T @ cell_velocity
        cofactor_rate = (
            np.trace(np.linalg.solve(cell, cell_velocity)) * geometry.cofactor
            - geometry.cofactor @ cell_velocity.T @ inverse_transpose
        )
        metric_rate = (
            cofactor_rate.T @ geometry.cofactor + geometry.cofactor.T @ cofactor_rate
        )
        geometric_force = inverse_transpose @ (
            np.trace(velocity_metric @ metric) * np.eye(3) - velocity_metric @ metric
        )
        acceleration = (
            driving + geometric_force - cell_velocity @ metric_rate
        ) @ inverse_metric
    elif model.mode.strain_cell:
        reference = cell_geometry(reference_cell)
        acceleration = driving @ (reference.metric / reference.volume**2)
    else:
        acceleration = driving
    # The original equations project strain acceleration for ALL variable-cell
    # modes. This is a continuous tensor operation, not space-group detection.
    strain_acceleration = np.linalg.solve(reference_cell.T, acceleration.T).T
    return 0.5 * (strain_acceleration + strain_acceleration.T) @ reference_cell
