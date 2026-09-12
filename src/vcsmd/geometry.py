"""Geometry operations for periodic cells.

The numerical core stores direct lattice vectors as the columns of a ``(3, 3)``
matrix.  Fractional coordinates are row vectors in the public API, so Cartesian
coordinates are obtained with ``fractional @ cell.T``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def wrap_fractional(positions: np.ndarray) -> np.ndarray:
    """Wrap into [0, 1), including roundoff from a tiny negative coordinate.

    IEEE remainder can round ``(-epsilon) % 1`` to exactly one. Canonicalizing
    that value to zero prevents an equivalent periodic point from violating
    the native state's half-open interval.
    """
    wrapped = np.remainder(positions, 1.0)
    return np.where(wrapped >= 1.0, 0.0, wrapped)


__all__ = [
    "CellGeometry",
    "cell_geometry",
    "cell_parameters",
    "replicate_cell",
    "validate_image_range",
]


def _readonly(array: np.ndarray) -> np.ndarray:
    """Return an owned, read-only float64 array."""

    result = np.array(array, dtype=np.float64, copy=True, order="C")
    result.setflags(write=False)
    return result


def _cell_array(cell: np.ndarray) -> np.ndarray:
    """Copy and validate the basic shape and finiteness of a cell."""

    try:
        result = np.array(cell, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "cell must be a finite numeric array with shape (3, 3)"
        ) from exc
    if result.shape != (3, 3):
        raise ValueError(f"cell must have shape (3, 3), got {result.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError("cell must contain only finite values")
    return result


def _validate_cell(cell: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """Validate a positively oriented, numerically nonsingular cell."""

    matrix = _cell_array(cell)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    largest = float(singular_values[0])
    smallest = float(singular_values[-1])
    # A relative criterion avoids accepting cells for which inversion is
    # dominated by roundoff while remaining independent of the unit scale.
    if (
        not np.isfinite(largest)
        or not np.isfinite(smallest)
        or smallest <= 1.0e-12 * largest
    ):
        raise ValueError("cell is singular or nearly singular")
    determinant = float(np.linalg.det(matrix))
    if not np.isfinite(determinant) or determinant <= 0.0:
        if determinant < 0.0:
            raise ValueError("cell has inverted (negative) orientation")
        raise ValueError("cell is singular or nearly singular")
    return matrix, determinant, singular_values


@dataclass(frozen=True, slots=True)
class CellGeometry:
    """Cached geometric tensors for a periodic cell.

    ``metric`` is ``cell.T @ cell`` and ``cofactor`` is the matrix satisfying
    ``cell.T @ cofactor = volume * I`` (equivalently
    ``volume * inv(cell).T``).  All array fields are independent read-only
    arrays.
    """

    volume: float
    metric: np.ndarray
    inverse_metric: np.ndarray
    cofactor: np.ndarray


def cell_geometry(cell: np.ndarray) -> CellGeometry:
    """Return volume, metric, inverse metric, and cofactor for ``cell``.

    The cell vectors are columns.  Inputs are copied and never modified.
    """

    matrix, volume, _ = _validate_cell(cell)
    metric = matrix.T @ matrix
    inverse_metric = np.linalg.inv(metric)
    cofactor = volume * np.linalg.inv(matrix).T
    return CellGeometry(
        volume=volume,
        metric=_readonly(metric),
        inverse_metric=_readonly(inverse_metric),
        cofactor=_readonly(cofactor),
    )


def cell_parameters(cell: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return vector lengths and ``(alpha, beta, gamma)`` angles in degrees.

    The angles are respectively between ``(b, c)``, ``(c, a)``, and ``(a, b)``.
    """

    matrix, _, _ = _validate_cell(cell)
    lengths = np.linalg.norm(matrix, axis=0)
    # alpha = angle(b,c), beta = angle(c,a), gamma = angle(a,b)
    pairs = ((1, 2), (2, 0), (0, 1))
    cosines = np.array(
        [
            np.dot(matrix[:, i], matrix[:, j]) / (lengths[i] * lengths[j])
            for i, j in pairs
        ],
        dtype=np.float64,
    )
    angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    return _readonly(lengths), _readonly(angles)


def _integer_triplet(
    values: object, *, name: str, strictly_positive: bool
) -> tuple[int, int, int]:
    """Validate an integer triplet without silently truncating floats."""

    try:
        array = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer sequence of length 3") from exc
    if array.shape != (3,) or array.dtype.kind not in "iu" or np.any(array < 0):
        qualifier = "positive" if strictly_positive else "nonnegative"
        raise ValueError(f"{name} must be a {qualifier} integer sequence of length 3")
    if strictly_positive and np.any(array == 0):
        raise ValueError(f"{name} must contain only positive integers")
    return tuple(int(x) for x in array)


def replicate_cell(
    cell: np.ndarray,
    fractional_positions: np.ndarray,
    repeats: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Replicate a cell and its basis in lexicographic shift order.

    For every shift ``(i, j, k)`` with ``i``, then ``j``, then ``k`` varying
    fastest in the nested loops, all original atoms are emitted in their input
    order.  The returned index map gives each emitted atom's source basis index.
    """

    matrix, _, _ = _validate_cell(cell)
    counts = _integer_triplet(repeats, name="repeats", strictly_positive=True)
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

    total = positions.shape[0] * counts[0] * counts[1] * counts[2]
    new_positions = np.empty((total, 3), dtype=np.float64)
    index_map = np.empty(total, dtype=np.intp)
    cursor = 0
    scale = np.asarray(counts, dtype=np.float64)
    for i in range(counts[0]):
        for j in range(counts[1]):
            for k in range(counts[2]):
                n = positions.shape[0]
                new_positions[cursor : cursor + n] = (positions + (i, j, k)) / scale
                index_map[cursor : cursor + n] = np.arange(n, dtype=np.intp)
                cursor += n

    new_cell = matrix @ np.diag(scale)
    new_cell.setflags(write=False)
    new_positions.setflags(write=False)
    index_map.setflags(write=False)
    return new_cell, new_positions, index_map


def validate_image_range(
    cell: np.ndarray,
    cutoff: float,
    image_range: tuple[int, int, int],
) -> None:
    """Ensure the requested image range reaches at least the cutoff radius.

    The required perpendicular face height for direction ``i`` is
    ``volume / area(face_i)``.  Thus ``n_i`` images in that direction cover a
    distance ``n_i * height_i``.  This check prevents silently truncating a
    periodic sum.
    """

    try:
        radius = float(cutoff)
    except (TypeError, ValueError) as exc:
        raise ValueError("cutoff must be a positive finite scalar") from exc
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("cutoff must be a positive finite scalar")
    counts = _integer_triplet(image_range, name="image_range", strictly_positive=False)
    matrix, _, _ = _validate_cell(cell)
    inverse = np.linalg.inv(matrix)
    heights = 1.0 / np.linalg.norm(inverse, axis=1)
    for direction, (count, height) in enumerate(zip(counts, heights)):
        covered = count * float(height)
        if covered < radius:
            raise ValueError(
                f"image_range direction {direction} is inadequate: "
                f"{count} * perpendicular face height {height:.16g} = {covered:.16g} "
                f"is less than cutoff {radius:.16g}"
            )
