# Original full workloads

The eight native workloads are complete original configurations. Each has TOML
and YAML forms; the first also has JSON. The matching historical input data are
linked in the table below. The source listings here show the first workload in
YAML, legacy format, and TOML.

To run the first workload directly from its native YAML configuration:

```bash
uv run vcsmd run examples/input-01.yaml --output-root runs
```

## Download the complete inputs

These downloads are copied directly from the repository during the documentation
build. They retain the original particle counts, physical parameters, and step
counts. The documentation build does not run or alter them. The legacy files
are Fortran-format input data for one-way import.

| Input | Mode | Particles | Steps | Download |
| --- | --- | ---: | ---: | --- |
| 1 | `metric-dyn` | 4 | 1,000 | {download}`TOML <../examples/input-01.toml>` · {download}`JSON <../examples/input-01.json>` · {download}`YAML <../examples/input-01.yaml>` · {download}`legacy <../examples/legacy/input-01.inp>` |
| 2 | `metric-dyn` | 4 | 2,000 | {download}`TOML <../examples/input-02.toml>` · {download}`YAML <../examples/input-02.yaml>` · {download}`legacy <../examples/legacy/input-02.inp>` |
| 3 | `metric-dyn` | 1 | 2,000 | {download}`TOML <../examples/input-03.toml>` · {download}`YAML <../examples/input-03.yaml>` · {download}`legacy <../examples/legacy/input-03.inp>` |
| 4 | `metric-dyn` | 1 | 2,000 | {download}`TOML <../examples/input-04.toml>` · {download}`YAML <../examples/input-04.yaml>` · {download}`legacy <../examples/legacy/input-04.inp>` |
| 5 | `metric-min` | 4 | 100 | {download}`TOML <../examples/input-05.toml>` · {download}`YAML <../examples/input-05.yaml>` · {download}`legacy <../examples/legacy/input-05.inp>` |
| 6 | `metric-min` | 2 | 100 | {download}`TOML <../examples/input-06.toml>` · {download}`YAML <../examples/input-06.yaml>` · {download}`legacy <../examples/legacy/input-06.inp>` |
| 7 | `fixed-dyn` | 32 | 1,000 | {download}`TOML <../examples/input-07.toml>` · {download}`YAML <../examples/input-07.yaml>` · {download}`legacy <../examples/legacy/input-07.inp>` |
| 8 | `fixed-dyn` | 32 | 1,000 | {download}`TOML <../examples/input-08.toml>` · {download}`YAML <../examples/input-08.yaml>` · {download}`legacy <../examples/legacy/input-08.inp>` |

The {download}`unmodified historical input collection <../examples/original-inputs.txt>`
is retained for provenance. The [legacy conversion guide](legacy-conversion.md)
explains the import boundary; the [validation report](validation.md) records
the full-run evidence and limitations.

## First workload in YAML

This native YAML file has the same configuration as the TOML and JSON forms.

```{literalinclude} ../examples/input-01.yaml
:language: yaml
```

## First workload in legacy format

This is the extracted Fortran-format input data for the first historical
workload. Use the [legacy conversion guide](legacy-conversion.md) to import it.

```{literalinclude} ../examples/legacy/input-01.inp
:language: text
```

## First workload in TOML

This listing is included from the source file so it stays aligned with the
download and the [getting-started commands](getting-started.md).

```{literalinclude} ../examples/input-01.toml
:language: toml
```
