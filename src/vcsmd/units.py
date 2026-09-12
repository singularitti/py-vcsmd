"""Quantities at the public boundary; numerical kernels do not import Pint."""

from __future__ import annotations

from typing import Any

import numpy as np
import pint

ureg = pint.UnitRegistry()
ureg.define("rydberg_time = hbar / rydberg")
ureg.define("rydberg_mass = rydberg * rydberg_time ** 2 / bohr ** 2")
Quantity = ureg.Quantity

LENGTH = "bohr"
ENERGY = "rydberg"
TIME = "rydberg_time"
MASS = "rydberg_mass"
PRESSURE = "rydberg / bohr ** 3"
VELOCITY = "bohr / rydberg_time"
FRACTIONAL_VELOCITY = "1 / rydberg_time"
CELL_INERTIA = "rydberg_mass / bohr ** 4"
BOLTZMANN = float(Quantity(1, "boltzmann_constant").to("rydberg / kelvin").magnitude)


def magnitude(value: Any, unit: str, *, name: str) -> np.ndarray:
    """Convert a quantity, including quantities from another registry.

    Dimensional bare numbers are deliberately rejected. SI base units provide
    an interchange for quantities from registries without our custom units.
    """
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{name} must be a Pint quantity with units")
    target = Quantity(1, unit)
    try:
        if value._REGISTRY is ureg:
            # Avoid a round trip through SI for our own quantities. Repeated
            # configuration serialization must not move normalized values by
            # a floating-point rounding unit on every read.
            local = value
        else:
            base = value.to_base_units()
            local = Quantity(base.magnitude, str(base.units))
        result = np.asarray(local.to(target.units).magnitude, dtype=np.float64)
    except (pint.PintError, ValueError, TypeError) as exc:
        raise ValueError(f"{name} must have units compatible with {unit}") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def scalar(value: Any, unit: str, *, name: str) -> float:
    result = magnitude(value, unit, name=name)
    if result.ndim != 0:
        raise ValueError(f"{name} must be a scalar quantity")
    return float(result)


def quantity_record(value: pint.Quantity) -> dict[str, Any]:
    """A format-independent, non-pickled quantity representation."""
    raw = np.asarray(value.magnitude)
    return {
        "value" if raw.ndim == 0 else "values": raw.tolist(),
        "unit": str(value.units),
    }


def parse_quantity(value: Any, *, name: str) -> pint.Quantity:
    if not isinstance(value, dict) or "unit" not in value:
        raise ValueError(f"{name} requires a value (or values) and a unit")
    keys = set(value)
    if keys not in ({"value", "unit"}, {"values", "unit"}):
        raise ValueError(f"{name} requires exactly {{value, unit}} or {{values, unit}}")
    try:
        raw = value.get("value", value.get("values"))
        return Quantity(np.asarray(raw, dtype=np.float64), value["unit"])
    except (pint.PintError, ValueError, TypeError) as exc:
        raise ValueError(f"Invalid quantity for {name}") from exc
