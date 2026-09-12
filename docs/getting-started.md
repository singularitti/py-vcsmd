# Getting started

VCSMD requires **Python 3.9 or newer**. The commands below assume a local
checkout and run from the repository root.

## Install

With [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run vcsmd --help
```

Alternatively, install the checkout into an existing Python environment:

```bash
python -m pip install .
python -m vcsmd --help
```

NumPy, Pint, PyYAML, and TOML support are installed as package dependencies.
The separate [documentation tools](documentation.md) require Python 3.12+.

## Choose an original configuration

The repository contains eight complete original workloads. The first is
available in equivalent TOML, JSON, and YAML forms; the
[examples page](examples.md) includes downloads and the full TOML source.

The first configuration contains four argon atoms, uses `metric-dyn`, and
requests all 1,000 original integration steps. Run it with:

```bash
uv run vcsmd run examples/input-01.toml --output-root runs
```

This command executes the full workload and writes a new run folder. Viewing
or building this documentation does not execute any of the displayed examples.

## Inspect the results

The application creates an hour-stamped, attempt-numbered folder under `runs/`.
Its README records the purpose, creation time, input provenance, and outcome.

```text
<run>/
  README.md
  inputs/                 configuration and copied implementation inputs
  outputs/                CSV data, metadata, and final checkpoint
    checkpoints/          periodic native checkpoints
```

CSV streams contain observables, cell history, particle trajectories, and
controller events. `outputs/checkpoint.npz` preserves the complete state needed
for native continuation. See [execution](reference/execution.md) for run options
and the [validation report](validation.md) for previously completed workloads.

## Work with Python objects

Use the functional interface when your application should control iteration
and output handling:

```python
from pathlib import Path

from vcsmd import initialize, prepare, simulate
from vcsmd.io import load_config

config = load_config(Path("examples/input-01.toml"))
model, initial_conditions = prepare(config)
state = initialize(model, initial_conditions, seed=config.seed)

for result in simulate(model, state, steps=config.steps):
    state = result.state
    temperature = result.observables.atomic_temperature.to("kelvin")
    print(result.observables.step_index, temperature)
```

`initialize` requires an explicit seed; native configurations default to 119.
The model and initial conditions contain normalized values. Public states and
observables expose physical quantities, and each step returns new, read-only
arrays. Iteration with `simulate` does not write files.

`result.observables` describes the completed integration step **before**
temperature rescaling or minimization quenching. `result.state` is the state
**after** those controller actions, ready for the next step or a checkpoint.

For application-managed outputs, use `vcsmd.execution.run` instead:

```python
from pathlib import Path

from vcsmd.execution import run
from vcsmd.io import load_config

source = Path("examples/input-01.toml")
report = run(
    load_config(source),
    output_root=Path("runs"),
    purpose="original argon workload",
    source_file=source,
)
print(report.status.value, report.completed_steps)
```

These are two ways to execute the same full workload. The
[public API reference](reference/public-api.md) documents the functional
interface; the [execution reference](reference/execution.md) covers file writing.

## Continue a native checkpoint

`--steps` is the number of **additional** steps:

```bash
uv run vcsmd resume runs/<run>/outputs/checkpoint.npz --steps 1000 --output-root runs
```

Replace `<run>` with the run-folder name. Native checkpoints preserve integrator
and controller history. Exact continuation is supported within the same
numerical software environment; different libraries or hardware can change
floating-point results. Historical partial states are not native checkpoints.

## Configuration conventions

JSON, YAML, and TOML share one schema:

- Dimensional scalars use `value` and `unit`; dimensional arrays use `values`
  and `unit`.
- Fractional positions are dimensionless arrays of shape `(N, 3)`. The cell is
  a dimensional `(3, 3)` matrix with cell vectors as its columns.
- The Lennard–Jones parameters are global for all particles. Use the `argon`
  preset or explicit `sigma` and `epsilon` records, never both.
- Variable-cell modes require `cell_inertia`. Ordinary cell modes use mass;
  metric and strain modes use mass divided by length to the fourth power.

| Formulation | Dynamics | Minimization |
| --- | --- | --- |
| Fixed cell | `fixed-dyn` | `fixed-min` |
| Parrinello–Rahman cell | `cell-dyn` | `cell-min` |
| Modified cell metric | `metric-dyn` | `metric-min` |
| Reference-cell strain | `strain-dyn` | `strain-min` |

Python callers use the corresponding `SimulationMode` enum members. See
[configuration](reference/configuration.md) for fields, defaults, enums, and
normalization, and [native I/O](reference/io.md) for serialization.
