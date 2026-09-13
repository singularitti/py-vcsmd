"""Periodic Lennard--Jones interactions for the numerical core."""

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np

from .geometry import validate_image_range, wrap_fractional

__all__ = ["ForceResult", "evaluate_lennard_jones"]


def _readonly(array: np.ndarray) -> np.ndarray:
    result = np.array(array, dtype=np.float64, copy=True, order="C")
    result.setflags(write=False)
    return result


@dataclass(frozen=True, slots=True)
class ForceResult:
    """Lennard--Jones energy, Cartesian forces, and configurational virial."""

    potential_energy: float
    forces: np.ndarray
    virial: np.ndarray


def _validate_inputs(
    cell: np.ndarray,
    fractional_positions: np.ndarray,
    sigma: float,
    epsilon: float,
    cutoff: float,
    image_range: tuple[int, int, int],
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray, float, float, float, tuple[int, int, int]]:
    try:
        matrix = np.array(cell, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError("cell must be a finite array with shape (3, 3)") from exc
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("cell must be a finite array with shape (3, 3)")
    try:
        positions = np.array(fractional_positions, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "fractional_positions must be a finite array with shape (N, 3)"
        ) from exc
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(
            f"fractional_positions must have shape (N, 3), got {positions.shape}"
        )
    if not np.all(np.isfinite(positions)):
        raise ValueError("fractional_positions must contain only finite values")
    # Fractional coordinates are periodic. Normalize the private copy once so
    # image validation and bounding-box pruning use one canonical box.
    positions = wrap_fractional(positions)
    try:
        sig = float(sigma)
        eps = float(epsilon)
        radius = float(cutoff)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "sigma, epsilon, and cutoff must be finite positive scalars"
        ) from exc
    if not all(np.isfinite(x) and x > 0.0 for x in (sig, eps, radius)):
        raise ValueError("sigma, epsilon, and cutoff must be finite positive scalars")
    if (
        isinstance(batch_size, (bool, np.bool_))
        or not isinstance(batch_size, (int, np.integer))
        or batch_size <= 0
    ):
        raise ValueError("batch_size must be a positive integer")
    # validate_image_range gives the detailed shape/type diagnostics and also
    # ensures no physically relevant image is omitted.
    validate_image_range(matrix, radius, image_range)
    counts = tuple(np.asarray(image_range).tolist())
    return matrix, positions, sig, eps, radius, tuple(int(x) for x in counts)


def _iter_valid_image_shift_batches(
    counts: tuple[int, int, int],
    position_minimum: np.ndarray,
    position_maximum: np.ndarray,
    singular_minimum: float,
    cutoff: float,
    output_batch_size: int,
) -> Iterator[np.ndarray]:
    """Yield retained translations in bounded vectorized batches."""

    n0, n1, n2 = (2 * count + 1 for count in counts)
    total = n0 * n1 * n2
    plane = n1 * n2
    pending: np.ndarray | None = None
    for start in range(0, total, output_batch_size):
        stop = min(start + output_batch_size, total)
        flat = np.arange(start, stop, dtype=np.int64)
        first = flat // plane - counts[0]
        remainder = flat % plane
        second = remainder // n2 - counts[1]
        third = remainder % n2 - counts[2]
        shifts = np.column_stack((first, second, third)).astype(np.float64, copy=False)
        # For delta = f_j - f_i + shift, each component lies in this interval.
        lower = position_minimum - position_maximum + shifts
        upper = position_maximum - position_minimum + shifts
        closest = np.where(lower > 0.0, lower, np.where(upper < 0.0, upper, 0.0))
        keep = singular_minimum * np.linalg.norm(closest, axis=1) <= cutoff
        # The caller handles the central-cell pairs separately. Excluding the
        # zero shift before buffering keeps noncentral batch boundaries exactly
        # aligned with the previous pair-major reduction order.
        keep &= np.any(shifts != 0.0, axis=1)
        retained = shifts[keep]
        if retained.size:
            if pending is None:
                pending = retained
            else:
                pending = np.concatenate((pending, retained), axis=0)
            while pending.shape[0] >= output_batch_size:
                yield pending[:output_batch_size]
                pending = pending[output_batch_size:]
                if pending.size == 0:
                    pending = None
                    break
    if pending is not None and pending.size:
        yield pending


def _pair_chunks(
    n_atoms: int,
    *,
    central_cell: bool,
    maximum_pairs: int,
) -> object:
    """Yield bounded ``(i, j)`` arrays without materializing all pairs."""

    if central_cell:
        # Keeping rows intact is useful for the common small-basis case while
        # ensuring that an arbitrarily large basis never creates N squared
        # temporary index arrays.
        for i in range(n_atoms):
            first_j = i + 1
            for j_start in range(first_j, n_atoms, maximum_pairs):
                j_stop = min(j_start + maximum_pairs, n_atoms)
                j = np.arange(j_start, j_stop, dtype=np.intp)
                yield np.full(j.shape, i, dtype=np.intp), j
        return

    rows_per_chunk = max(1, maximum_pairs // max(1, n_atoms))
    if n_atoms > maximum_pairs:
        for i in range(n_atoms):
            for j_start in range(0, n_atoms, maximum_pairs):
                j_stop = min(j_start + maximum_pairs, n_atoms)
                j = np.arange(j_start, j_stop, dtype=np.intp)
                yield np.full(j.shape, i, dtype=np.intp), j
        return
    for i_start in range(0, n_atoms, rows_per_chunk):
        i_stop = min(i_start + rows_per_chunk, n_atoms)
        i = np.arange(i_start, i_stop, dtype=np.intp)
        j = np.arange(n_atoms, dtype=np.intp)
        yield np.repeat(i, n_atoms), np.tile(j, i.size)


def evaluate_lennard_jones(
    cell: np.ndarray,
    fractional_positions: np.ndarray,
    *,
    sigma: float,
    epsilon: float,
    cutoff: float,
    image_range: tuple[int, int, int],
    batch_size: int = 65536,
) -> ForceResult:
    """Evaluate an unshifted, inclusive-cutoff Lennard--Jones potential.

    The central cell's distinct pairs are counted once.  Every nonzero image
    shift is evaluated for all ordered particle pairs; those terms receive a
    half weight for energy and virial, while each central particle receives
    the full force from every image.  This retains self-image interactions and
    gives the periodic force and virial their usual translationally invariant
    values without constructing an all-pairs/all-images tensor.
    """

    matrix, positions, sig, eps, radius, counts = _validate_inputs(
        cell,
        fractional_positions,
        sigma,
        epsilon,
        cutoff,
        image_range,
        batch_size,
    )
    n_atoms = positions.shape[0]
    forces = np.zeros((n_atoms, 3), dtype=np.float64)
    virial = np.zeros((3, 3), dtype=np.float64)
    if n_atoms == 0:
        return ForceResult(0.0, _readonly(forces), _readonly(virial))

    singular_minimum = float(np.linalg.svd(matrix, compute_uv=False)[-1])
    position_minimum = positions.min(axis=0)
    position_maximum = positions.max(axis=0)
    cutoff_squared = radius * radius
    sigma_squared = sig * sig
    potential_energy = 0.0

    def accumulate_pair_batches(
        pair_i: np.ndarray,
        pair_j: np.ndarray,
        image_shifts: np.ndarray,
        *,
        central_cell: bool,
    ) -> None:
        nonlocal potential_energy, virial
        n_pairs = pair_i.size
        if n_pairs == 0 or image_shifts.size == 0:
            return
        # Combining pair chunks with this many images keeps every flattened
        # image/pair batch bounded by batch_size.
        images_per_batch = max(1, batch_size // n_pairs)
        for image_start in range(0, image_shifts.shape[0], images_per_batch):
            image_batch = image_shifts[image_start : image_start + images_per_batch]
            n_images = image_batch.shape[0]
            flat_i = np.tile(pair_i, n_images)
            flat_j = np.tile(pair_j, n_images)
            flat_shifts = np.repeat(image_batch, n_pairs, axis=0)
            delta_fractional = positions[flat_j] - positions[flat_i] + flat_shifts
            displacement = delta_fractional @ matrix.T
            distance_squared = np.einsum(
                "ij,ij->i", displacement, displacement, optimize=True
            )
            selected = distance_squared <= cutoff_squared
            if not np.any(selected):
                continue
            selected_indices = np.flatnonzero(selected)
            selected_distance_squared = distance_squared[selected_indices]
            if np.any(selected_distance_squared <= 0.0):
                collision_kind = (
                    "distinct central particles"
                    if central_cell
                    else "central particle and periodic image"
                )
                raise ValueError(
                    f"zero-distance Lennard-Jones collision between {collision_kind}"
                )
            selected_displacement = displacement[selected_indices]
            inv_r2 = sigma_squared / selected_distance_squared
            inv_r6 = inv_r2 * inv_r2 * inv_r2
            inv_r12 = inv_r6 * inv_r6
            pair_energy = 4.0 * eps * (inv_r12 - inv_r6)
            factor = 24.0 * eps / selected_distance_squared * (2.0 * inv_r12 - inv_r6)
            pair_force = factor[:, None] * selected_displacement
            pair_weight = 1.0 if central_cell else 0.5
            potential_energy += pair_weight * float(np.sum(pair_energy))
            virial += pair_weight * np.einsum(
                "ni,nj->ij", pair_force, selected_displacement, optimize=True
            )
            selected_i = flat_i[selected_indices]
            for axis in range(3):
                forces[:, axis] -= np.bincount(
                    selected_i,
                    weights=pair_force[:, axis],
                    minlength=n_atoms,
                )
            if central_cell:
                selected_j = flat_j[selected_indices]
                for axis in range(3):
                    forces[:, axis] += np.bincount(
                        selected_j,
                        weights=pair_force[:, axis],
                        minlength=n_atoms,
                    )

    central_shifts = np.zeros((1, 3), dtype=np.float64)
    for pair_i, pair_j in _pair_chunks(
        n_atoms,
        central_cell=True,
        maximum_pairs=batch_size,
    ):
        accumulate_pair_batches(pair_i, pair_j, central_shifts, central_cell=True)
    for pair_i, pair_j in _pair_chunks(
        n_atoms,
        central_cell=False,
        maximum_pairs=batch_size,
    ):
        images_per_batch = max(1, batch_size // pair_i.size)
        for noncentral_shifts in _iter_valid_image_shift_batches(
            counts,
            position_minimum,
            position_maximum,
            singular_minimum,
            radius,
            images_per_batch,
        ):
            if noncentral_shifts.size:
                accumulate_pair_batches(
                    pair_i,
                    pair_j,
                    noncentral_shifts,
                    central_cell=False,
                )

    return ForceResult(
        potential_energy=float(potential_energy),
        forces=_readonly(forces),
        virial=_readonly(virial),
    )
