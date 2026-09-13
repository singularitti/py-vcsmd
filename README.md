# Python VCSMD

`vcsmd` is a Python 3.10+ implementation of classical variable-cell molecular dynamics with a pure NumPy numerical core and a Pint-aware application boundary. It implements the eight active fixed-cell, Parrinello–Rahman, modified-metric, and reference-strain formulations described by the historical VCSMD project. The numerical model contains resolved arrays and scalars; configuration files, checkpoints, legacy files, progress, and output writing live in adapters around that core.

The scientific scope is hydrostatic-pressure classical Lennard–Jones dynamics. The active potential is configurable through global `sigma` and `epsilon` parameters, with the original argon values available as the named `argon` preset. Tabulated potentials, electronic-structure calculations, tensor external stress, crystal-standardization operations, and inactive historical experiments are outside this package.

To understand the implementation, start with the [code structure and walkthrough](docs/architecture.md). It explains the objects, module responsibilities, units, and the sequence of one integration step. The [scientific guide](docs/scientific-guide.md) develops the equations, and the [validation report](docs/validation.md) records what has been checked.

## Install

The package requires Python 3.10 or newer. With [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run vcsmd --help
```

The package can also be installed into an existing compatible environment with `pip install .`. NumPy, Pint, PyYAML, and TOML support are installed from the project dependencies.

## Documentation website

The documentation site uses Sphinx, MyST Markdown, and Furo, with searchable
Python API references, rendered equations, and architecture diagrams. Building
the docs requires Python 3.12 or newer:

```bash
uv sync --group docs --no-default-groups --python 3.14
uv run --group docs --no-default-groups sphinx-build -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`, or use the live preview command in the
[documentation build guide](docs/documentation.md). The build reads docstrings
and the original example files without executing simulations.

## Simulation modes

Python code uses `SimulationMode` members. Native configuration files and command-line workflows use the serialized values below.

| Scientific formulation | Python member | Serialized value |
| --- | --- | --- |
| Fixed-cell dynamics | `SimulationMode.FIXED_DYNAMICS` | `fixed-dyn` |
| Fixed-cell minimization | `SimulationMode.FIXED_MINIMIZATION` | `fixed-min` |
| Parrinello–Rahman cell dynamics | `SimulationMode.CELL_DYNAMICS` | `cell-dyn` |
| Parrinello–Rahman cell minimization | `SimulationMode.CELL_MINIMIZATION` | `cell-min` |
| Modified cell-metric dynamics | `SimulationMode.METRIC_DYNAMICS` | `metric-dyn` |
| Modified cell-metric minimization | `SimulationMode.METRIC_MINIMIZATION` | `metric-min` |
| Reference-cell strain dynamics | `SimulationMode.STRAIN_DYNAMICS` | `strain-dyn` |
| Reference-cell strain minimization | `SimulationMode.STRAIN_MINIMIZATION` | `strain-min` |

Initialization behavior is also an enum: `InitializationMode.THERMAL` or `InitializationMode.PROVIDED`.

## Python interface

Use the top-level `vcsmd` API to parse a native file, resolve quantities and formulation-specific units, and evolve an immutable state:

```python
from pathlib import Path

from vcsmd import initialize, load_config, prepare, simulate

config = load_config(Path("examples/input-01.toml"))
model, initial_conditions = prepare(config)
state = initialize(model, initial_conditions, seed=config.seed)

for result in simulate(model, state, steps=config.steps):
    state = result.state
    temperature = result.observables.atomic_temperature
    print(result.observables.step_index, temperature)
```

The normalized `NumericalModel` and `InitialConditions` are immutable and use float64 arrays in bohr, Rydberg, Rydberg time, and the corresponding mass units. The public `SimulationState` and `Observables` returned by `vcsmd.initialize`, `vcsmd.step`, and `vcsmd.simulate` expose dimensional fields as Pint quantities. Published arrays are read-only and each evolution step owns its returned arrays.

For an interactive session, the package also supports wildcard imports:

```python
from pathlib import Path

from vcsmd import *

config = load_config(Path("examples/input-01.toml"))
model, initial_conditions = prepare(config)
state = initialize(model, initial_conditions, seed=config.seed)
for result in simulate(model, state, steps=config.steps):
    state = result.state
```

`vcsmd.__all__` defines the names included by `from vcsmd import *`. Explicit
imports are recommended in maintained application and library code because they
make dependencies visible; the star form is supported as a convenient interactive
workflow.

For a complete application run, use the file-writing orchestration layer:

```python
from pathlib import Path

from vcsmd import load_config, run

report = run(
    load_config(Path("examples/input-01.toml")),
    output_root=Path("runs"),
    purpose="argon example",
)
print(report.status.value, report.directory)
```

The numerical functions do not access files, print, log, inspect the environment, or mutate global random state. The application layer owns those concerns.

## Native configuration

JSON, YAML, and TOML use one schema and the same validation pipeline. Dimensional scalars are records with `value` and `unit`; dimensional arrays use `values` and `unit`. The first original full workload is provided in equivalent [TOML](examples/input-01.toml), [JSON](examples/input-01.json), and [YAML](examples/input-01.yaml) configurations. It contains four atoms, uses `metric-dyn`, and runs all 1,000 original steps. All eight unchanged workloads and their provenance are described in the [examples guide](examples/README.md).

The potential record accepts either `preset = "argon"` or explicit dimensional `sigma` and `epsilon` records. These alternatives cannot be combined. Fractional positions are dimensionless `(N, 3)` arrays; the cell is a dimensional `(3, 3)` matrix whose columns are the cell vectors. The application default seed is 119.

Variable-cell modes require `cell_inertia`. Ordinary cell modes use mass units; modified-metric and reference-strain modes use mass/length⁴. A provided initialization adds `fractional_velocities` and may add `cell_velocity`, each as dimensional arrays. Provided particle velocities must have zero mass-weighted global center-of-mass momentum. `save_config` and `load_config` select JSON, YAML, or TOML from the `pathlib.Path` suffix; YAML is loaded safely.

## Command line

Run a native configuration:

```bash
uv run vcsmd run examples/input-01.toml --output-root runs
```

Resume a complete native checkpoint for a specified number of additional steps:

```bash
uv run vcsmd resume runs/<run>/outputs/checkpoint.npz --steps 1000 --output-root runs
```

Import one historical input into a native configuration, or convert available historical output files into native datasets and a conversion report:

```bash
uv run vcsmd import-legacy <legacy-input> examples/imported.json
uv run vcsmd convert-legacy <legacy-run-directory>
```

Legacy conversion is one-way. It recognizes historical filenames, fixed-width records, Fortran exponents, historical unit conventions, mode codes, and coordinate-operation flags only inside `vcsmd.compat`. By default, the command creates an hour-stamped, attempt-numbered conversion directory under `runs/`. An explicit `--destination` must be empty and outside the source directory. Imported CSV datasets preserve the fields actually available in historical files; their schema is documented in the [compatibility guide](docs/legacy-conversion.md). Incomplete historical state is reported as incomplete data and is never presented as an exact native checkpoint.

## Run folders and outputs

The top-level `vcsmd.run` function creates an hour-stamped, attempt-numbered directory containing `inputs/`, `outputs/`, and `outputs/checkpoints/`. The folder README records its purpose, creation time, inputs, provenance, outputs, and completion status. `run_metadata.json` records the schema, code and dependency hashes, configuration provenance, units, requested and completed steps, and failure information without exposing machine-specific paths.

Native outputs include unit-labelled CSV streams for observables, cell history, trajectories, and controller events. `checkpoint.npz` and periodic checkpoint archives use named non-pickled arrays and preserve the cell, positions, velocities, current and previous accelerations, reference cell, running accumulators, controller counter, and all model parameters needed for continuation. Observables describe the completed integration step before temperature rescaling and minimization quenching; trajectories and checkpoints represent the state after those events.

Native checkpoints preserve exact continuation within the same numerical software environment. Floating-point libraries, hardware, or package changes may prevent bitwise equality across environments.

## Verification

All eight supplied workloads completed at their original sizes, parameters, and step counts: 9,200 integration steps. Continuing the second workload from its step-1,000 checkpoint to step 2,000 produced exactly the same complete final state as the uninterrupted run. Post-feature checks covered force and energy derivatives, periodic virial, cell-metric equations, deterministic initialization, input immutability, units, and format equivalence. Installation and imports were verified on Python 3.10 and 3.14.

The fixtures exercise `fixed-dyn`, `metric-dyn`, and `metric-min`. The five other modes are implemented but have no suitable supplied full-trajectory cases. An executed Fortran reference and complete historical output fixtures were unavailable. Completion and finite values do not establish convergence or energy conservation for the original parameters.

The [scientific guide](docs/scientific-guide.md) documents equations, unit conventions, numerical boundaries, and corrected historical defects. The [validation report](docs/validation.md) records results, preserved CSV/checkpoint/PDF artifacts, and the limits of the evidence.

The [release guide](docs/releasing.md) explains the CI checks, PyPI Trusted Publishing setup, and package version updates required for publication.
