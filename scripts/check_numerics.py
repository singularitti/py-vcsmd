"""Post-feature numerical identities using complete existing workload data.

No trajectory is integrated here and no synthetic or reduced simulation is
created. Derivative probes evaluate the full periodic system from saved runs.
"""

from __future__ import annotations

import argparse
import ast
import csv
from dataclasses import replace
from pathlib import Path

import numpy as np

import vcsmd
from vcsmd.cell_dynamics import cell_acceleration
from vcsmd.config import prepare
from vcsmd.dynamics import initialize
from vcsmd.geometry import cell_geometry, wrap_fractional
from vcsmd.io import load_config
from vcsmd.io.checkpoint import load_checkpoint
from vcsmd.models import RunningAverages
from vcsmd.potentials import evaluate_lennard_jones


def _force(model, cell, positions):
    return evaluate_lennard_jones(
        cell,
        positions,
        sigma=model.sigma,
        epsilon=model.epsilon,
        cutoff=model.cutoff,
        image_range=model.image_range,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("validation_run", type=Path)
    args = parser.parse_args()
    output = args.validation_run / "outputs"
    with (output / "summary.csv").open(newline="") as stream:
        runs = {
            int(row["input"]): args.validation_run / row["run"]
            for row in csv.DictReader(stream)
        }
    rows = []

    def record(name: str, passed: bool, error: float = 0.0) -> None:
        rows.append(
            {"check": name, "passed": bool(passed), "maximum_error": float(error)}
        )
        print(name, "PASS" if passed else "FAIL", f"({error:.3g})", flush=True)

    config = load_config(args.validation_run / "inputs" / "input-07.json")
    model, initial = prepare(config)
    source_arrays = (
        initial.cell.copy(),
        initial.fractional_positions.copy(),
        model.masses.copy(),
    )
    global_before = np.random.get_state()
    first = initialize(model, initial, seed=config.seed)
    second = initialize(model, initial, seed=config.seed)
    global_after = np.random.get_state()
    record(
        "Seeded full-system initialization is deterministic",
        np.array_equal(first.fractional_velocities, second.fractional_velocities),
    )
    record(
        "Initialization preserves global random state",
        all(
            np.array_equal(a, b) if isinstance(a, np.ndarray) else a == b
            for a, b in zip(global_before, global_after)
        ),
    )
    record(
        "Initialization leaves supplied arrays unchanged",
        all(
            np.array_equal(a, b)
            for a, b in zip(
                source_arrays,
                (initial.cell, initial.fractional_positions, model.masses),
            )
        ),
    )
    record(
        "Published numerical arrays are read-only",
        not first.cell.flags.writeable
        and not first.fractional_positions.flags.writeable
        and not model.masses.flags.writeable,
    )
    velocity = first.fractional_velocities @ first.cell.T
    momentum = np.sum(model.masses[:, None] * velocity, axis=0)
    record(
        "Thermal initialization removes global momentum",
        np.max(np.abs(momentum)) < 1e-9,
        np.max(np.abs(momentum)),
    )
    kinetic = 0.5 * np.sum(model.masses[:, None] * velocity**2)
    temperature = 2 * kinetic / (3 * (len(velocity) - 1) * model.boltzmann_constant)
    record(
        "Thermal initialization reaches requested temperature",
        abs(temperature - model.temperature) < 1e-9,
        abs(temperature - model.temperature),
    )

    model, state = load_checkpoint(
        runs[7] / "outputs" / "checkpoints" / "step-00000500.npz"
    )
    cell = state.cell
    positions = state.fractional_positions
    analytical = _force(model, cell, positions)
    derivative = np.empty_like(positions)
    h = 1e-5  # Cartesian bohr; evaluate every component of the FULL saved system.
    inverse = np.linalg.inv(cell)
    for atom in range(len(positions)):
        for axis in range(3):
            plus = positions.copy()
            minus = positions.copy()
            plus[atom] += h * inverse[:, axis]
            minus[atom] -= h * inverse[:, axis]
            derivative[atom, axis] = -(
                _force(model, cell, plus).potential_energy
                - _force(model, cell, minus).potential_energy
            ) / (2 * h)
    error = float(np.max(np.abs(derivative - analytical.forces)))
    record(
        "Full-system force equals negative energy gradient",
        np.allclose(derivative, analytical.forces, atol=1e-9, rtol=1e-5),
        error,
    )
    record(
        "Periodic forces conserve total momentum",
        np.max(np.abs(np.sum(analytical.forces, axis=0))) < 1e-11,
        np.max(np.abs(np.sum(analytical.forces, axis=0))),
    )

    model, state = load_checkpoint(
        runs[3] / "outputs" / "checkpoints" / "step-00001000.npz"
    )
    cell, rate = state.cell, state.cell_velocity
    forces = _force(model, cell, state.fractional_positions)
    geometry = cell_geometry(cell)
    record(
        "Cofactor orientation and inverse metric",
        np.allclose(cell.T @ geometry.cofactor, geometry.volume * np.eye(3), atol=1e-10)
        and np.allclose(
            geometry.metric @ geometry.inverse_metric, np.eye(3), atol=1e-12
        ),
    )
    # The virtual deformations do not integrate a new configuration: they
    # differentiate the saved state's full periodic potential and kinetic form.
    deformation_gradient = np.empty((3, 3))
    dh = 1e-6
    for i in range(3):
        for j in range(3):
            deformation = np.zeros((3, 3))
            deformation[i, j] = dh
            positive = _force(
                model, (np.eye(3) + deformation) @ cell, state.fractional_positions
            )
            negative = _force(
                model, (np.eye(3) - deformation) @ cell, state.fractional_positions
            )
            deformation_gradient[i, j] = (
                positive.potential_energy - negative.potential_energy
            ) / (2 * dh)
    error = float(np.max(np.abs(deformation_gradient + forces.virial)))
    record(
        "Periodic virial equals negative deformation gradient",
        np.allclose(deformation_gradient, -forces.virial, atol=1e-8, rtol=1e-5),
        error,
    )

    def metric(matrix):
        cofactor = cell_geometry(matrix).cofactor
        return cofactor.T @ cofactor

    def kinetic_form(matrix):
        return 0.5 * float(np.trace(rate @ metric(matrix) @ rate.T))

    kinetic_gradient = np.empty((3, 3))
    delta = 1e-5
    for i in range(3):
        for j in range(3):
            perturbation = np.zeros((3, 3))
            perturbation[i, j] = delta
            kinetic_gradient[i, j] = (
                kinetic_form(cell + perturbation) - kinetic_form(cell - perturbation)
            ) / (2 * delta)
    time_delta = 1e-5 * np.linalg.norm(cell) / max(np.linalg.norm(rate), 1e-15)
    metric_rate = (
        metric(cell + time_delta * rate) - metric(cell - time_delta * rate)
    ) / (2 * time_delta)
    velocity = state.fractional_velocities @ cell.T
    stress = ((velocity.T * model.masses) @ velocity + forces.virial) / geometry.volume
    driving = (
        (stress - model.external_pressure * np.eye(3))
        @ geometry.cofactor
        / model.cell_inertia
    )
    unprojected = (driving + kinetic_gradient - rate @ metric_rate) @ np.linalg.inv(
        metric(cell)
    )
    strain = unprojected @ np.linalg.inv(state.reference_cell)
    expected = 0.5 * (strain + strain.T) @ state.reference_cell
    computed = cell_acceleration(
        model, cell, rate, state.reference_cell, velocity, forces.virial
    )
    error = float(np.max(np.abs(expected - computed)))
    record(
        "Metric acceleration agrees with differentiated Lagrangian",
        np.allclose(expected, computed, atol=1e-11, rtol=1e-5),
        error,
    )

    metric_config = load_config(args.validation_run / "inputs" / "input-01.json")
    rejected = False
    try:
        prepare(replace(metric_config, cell_inertia=vcsmd.Quantity(1, "kilogram")))
    except ValueError:
        rejected = True
    record("Metric inertia rejects ordinary mass dimensions", rejected)
    record(
        "Periodic wrapping respects half-open interval",
        np.all(
            wrap_fractional(np.nextafter(state.fractional_positions, -np.inf)) < 1.0
        ),
    )
    for invalid in (-1, 0.5, True):
        rejected = False
        try:
            RunningAverages(count=invalid)
        except ValueError:
            rejected = True
        record(f"Reject invalid running-average count {invalid!r}", rejected)
    package = Path(vcsmd.__file__).parent
    forbidden = {
        "os",
        "pathlib",
        "pint",
        "io",
        "compat",
        "execution",
        "cli",
        "json",
        "yaml",
    }
    violations = []
    for module in (
        "models",
        "geometry",
        "potentials",
        "cell_dynamics",
        "dynamics",
        "observables",
    ):
        tree = ast.parse((package / f"{module}.py").read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports = [item.name for item in node.names]
            elif isinstance(node, ast.ImportFrom):
                imports = [node.module or ""]
            else:
                continue
            if any(name.split(".")[0] in forbidden for name in imports):
                violations.append(module)
    record("Numerical modules do not import I/O, legacy or Pint", not violations)
    with (output / "numerical_checks.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["check", "passed", "maximum_error"])
        writer.writeheader()
        writer.writerows(rows)
    return 0 if all(row["passed"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
