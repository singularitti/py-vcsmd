"""Unit-aware, functional variable-cell molecular dynamics.

Import the native API directly from this package. ``__all__`` defines the
supported wildcard-import namespace; explicit imports are clearer in maintained
code. Implementations import their dependencies from concrete submodules rather
than importing objects back from this facade.
"""

from .config import LennardJones, SimulationConfig, Structure, prepare
from .execution import RunReport, RunStatus, resume, run
from .io import (
    config_from_mapping,
    config_to_mapping,
    load_checkpoint,
    load_config,
    save_checkpoint,
    save_config,
)
from .models import (
    EventKind,
    InitialConditions,
    InitializationMode,
    NumericalModel,
    SimulationEvent,
    SimulationMode,
    TemperatureControl,
)
from .simulation import (
    Observables,
    SimulationState,
    StepResult,
    initialize,
    simulate,
    step,
    with_units,
)
from .units import Quantity, ureg

__all__ = [
    "EventKind",
    "InitialConditions",
    "InitializationMode",
    "LennardJones",
    "NumericalModel",
    "Observables",
    "Quantity",
    "RunReport",
    "RunStatus",
    "SimulationConfig",
    "SimulationEvent",
    "SimulationMode",
    "SimulationState",
    "StepResult",
    "Structure",
    "TemperatureControl",
    "config_from_mapping",
    "config_to_mapping",
    "initialize",
    "load_checkpoint",
    "load_config",
    "prepare",
    "resume",
    "run",
    "save_checkpoint",
    "save_config",
    "simulate",
    "step",
    "ureg",
    "with_units",
]

__version__ = "0.1.0"
