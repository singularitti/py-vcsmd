# Original full workloads

These are native conversions of the eight complete examples in the original input collection. Atom counts, step counts and numerical parameters are preserved; no reduced or synthetic simulation is included.

The first workload is provided in JSON, YAML and TOML to demonstrate equivalent serialization. All use the active Lennard–Jones potential. The unused table appended to the historical eighth example is retained only in the unmodified source collection.

Input provenance: copied from the original VCSMD repository input collection. SHA-256: f7178e7004d731883f305f66945c568d3e5129aacd00fa611e8bdf8ae52b91f1.

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
the same functions, and the configuration parameters are unchanged.
