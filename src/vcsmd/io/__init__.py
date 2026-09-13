"""Native configuration and checkpoint file adapters.

The numerical package has no knowledge of file formats.  This namespace is the
single boundary for native JSON, YAML, TOML, and checkpoint files.
"""

from pathlib import Path

from ..models import NumericalModel
from ..simulation import SimulationState
from . import checkpoint
from .config import config_from_mapping, config_to_mapping, load_config, save_config


def load_checkpoint(path: Path) -> tuple[NumericalModel, SimulationState]:
    """Load a model and a quantity-bearing public state for ``vcsmd.simulate``."""
    model, state = checkpoint.load_checkpoint(path)
    return model, SimulationState(state)


def save_checkpoint(path: Path, model: NumericalModel, state: SimulationState) -> None:
    """Persist a public simulation state, including its complete numerical history."""
    if not isinstance(state, SimulationState):
        raise TypeError(
            "state must be the public SimulationState returned by initialize or step"
        )
    checkpoint.save_checkpoint(path, model, state.numerical)


__all__ = [
    "config_from_mapping",
    "config_to_mapping",
    "load_checkpoint",
    "load_config",
    "save_checkpoint",
    "save_config",
]
