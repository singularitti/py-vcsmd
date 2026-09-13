# Original full workloads

These are native conversions of the eight complete examples in the original input collection. Atom counts, step counts and numerical parameters are preserved; no reduced or synthetic simulation is included. Every workload is available as TOML and YAML; the first is also available as JSON. The matching historical Fortran-format input data are in [`legacy/`](legacy/), including [workload 1](legacy/input-01.inp). The [examples guide](../docs/examples.md) has the complete download table and first-workload source listings.

The first workload is provided in JSON, YAML and TOML to demonstrate equivalent serialization. All use the active Lennard–Jones potential. The unused table appended to the historical eighth example is retained only in the unmodified source collection. It is not part of the individual extracted legacy input files.

Input provenance: copied from the original VCSMD repository input collection. SHA-256: f7178e7004d731883f305f66945c568d3e5129aacd00fa611e8bdf8ae52b91f1.

## Download the inputs

| Workload | TOML | YAML | Legacy input data |
| --- | --- | --- | --- |
| 1 | [input-01.toml](input-01.toml) | [input-01.yaml](input-01.yaml) | [input-01.inp](legacy/input-01.inp) |
| 2 | [input-02.toml](input-02.toml) | [input-02.yaml](input-02.yaml) | [input-02.inp](legacy/input-02.inp) |
| 3 | [input-03.toml](input-03.toml) | [input-03.yaml](input-03.yaml) | [input-03.inp](legacy/input-03.inp) |
| 4 | [input-04.toml](input-04.toml) | [input-04.yaml](input-04.yaml) | [input-04.inp](legacy/input-04.inp) |
| 5 | [input-05.toml](input-05.toml) | [input-05.yaml](input-05.yaml) | [input-05.inp](legacy/input-05.inp) |
| 6 | [input-06.toml](input-06.toml) | [input-06.yaml](input-06.yaml) | [input-06.inp](legacy/input-06.inp) |
| 7 | [input-07.toml](input-07.toml) | [input-07.yaml](input-07.yaml) | [input-07.inp](legacy/input-07.inp) |
| 8 | [input-08.toml](input-08.toml) | [input-08.yaml](input-08.yaml) | [input-08.inp](legacy/input-08.inp) |

The [validation report](../docs/validation.md) records the format-coverage
run: all eight YAML and imported-legacy configurations matched their TOML
counterparts in configuration comparison, and each YAML/legacy pair matched in
final checkpoints and complete CSV outputs.

The first workload is also available as [JSON](input-01.json). To import its
legacy input and run the resulting native configuration:

```bash
mkdir -p runs
uv run vcsmd import-legacy examples/legacy/input-01.inp runs/imported-input-01.yaml
uv run vcsmd run runs/imported-input-01.yaml --output-root runs
```

## Run a full workload from Python

From the repository root, the package-root API can run the complete first input:

```python
from pathlib import Path

from vcsmd import *

source = Path("examples/input-01.toml")
report = run(
    load_config(source),
    source_file=source,
    output_root=Path("runs"),
    purpose="original-input-1",
)
```

This executes all 1,000 original steps and writes native results. The wildcard
form is convenient in an interactive session. In maintained code, use the
explicit equivalent `from vcsmd import load_config, run`. Both forms expose
the same functions, and the configuration parameters are unchanged. To import
the corresponding historical input data first, see the [legacy conversion
guide](../docs/legacy-conversion.md).
