```{include} ../examples/README.md
```

## Download the complete inputs

These downloads are copied directly from the repository during the documentation
build. They retain the original particle counts, physical parameters, and step
counts. The documentation build does not run or alter them.

| Input | Mode | Particles | Steps | Download |
| --- | --- | ---: | ---: | --- |
| 1 | `metric-dyn` | 4 | 1,000 | {download}`TOML <../examples/input-01.toml>` · {download}`JSON <../examples/input-01.json>` · {download}`YAML <../examples/input-01.yaml>` |
| 2 | `metric-dyn` | 4 | 2,000 | {download}`TOML <../examples/input-02.toml>` |
| 3 | `metric-dyn` | 1 | 2,000 | {download}`TOML <../examples/input-03.toml>` |
| 4 | `metric-dyn` | 1 | 2,000 | {download}`TOML <../examples/input-04.toml>` |
| 5 | `metric-min` | 4 | 100 | {download}`TOML <../examples/input-05.toml>` |
| 6 | `metric-min` | 2 | 100 | {download}`TOML <../examples/input-06.toml>` |
| 7 | `fixed-dyn` | 32 | 1,000 | {download}`TOML <../examples/input-07.toml>` |
| 8 | `fixed-dyn` | 32 | 1,000 | {download}`TOML <../examples/input-08.toml>` |

The {download}`unmodified historical input collection <../examples/original-inputs.txt>`
is retained for provenance. The [legacy conversion guide](legacy-conversion.md)
explains the import boundary; the [validation report](validation.md) records
the full-run evidence and limitations.

## First workload in TOML

This listing is included from the source file so it stays aligned with the
download and the [getting-started commands](getting-started.md).

```{literalinclude} ../examples/input-01.toml
:language: toml
```
