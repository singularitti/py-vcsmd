# VCSMD

**Classical variable-cell molecular dynamics in Python.**

VCSMD combines a pure NumPy numerical core with a unit-aware Python interface.
Describe a periodic structure with Pint quantities, choose a fixed-cell or
variable-cell formulation, and evolve immutable states or write a complete,
documented run to disk.

The package implements eight dynamics and minimization modes with a global
Lennard–Jones potential and hydrostatic external pressure. It supports JSON,
YAML, and TOML configurations, CSV results, native checkpoints, and one-way
conversion of historical VCSMD data.

## Start here

| What you want to do | Where to start |
| --- | --- |
| Install VCSMD and run an original workload | [Getting started](getting-started.md) |
| Download the complete original inputs | [Example configurations](examples.md) |
| Understand units, state, and the integration step | [Implementation walkthrough](architecture.md) |
| Read the equations and physical assumptions | [Scientific guide](scientific-guide.md) |
| Look up Python functions and data models | [API reference](reference/index.md) |
| Assess the available numerical evidence | [Validation report](validation.md) |

## Functional core, explicit boundaries

The Python workflow is `load_config → prepare → initialize → simulate`.
Configuration objects retain physical units. Normalized models supply the
numerical solver with arrays and scalars. Public results expose dimensional
fields as Pint quantities, and published arrays are read-only.

The execution layer owns file writing, checkpoints, progress, and provenance.
Each application run records its inputs and outputs in a dedicated folder.
See the [execution reference](reference/execution.md) for that boundary.

:::{note}
The supplied full workloads exercise `fixed-dyn`, `metric-dyn`, and `metric-min`.
The other five modes are implemented but lack supplied full-trajectory cases.
Electronic-structure calculations, tensor external stress, tabulated potentials,
and crystal standardization are outside the implemented scope. The
[validation report](validation.md) describes the evidence and its limits.
:::

```{toctree}
:hidden:
:caption: User guide

getting-started
examples
legacy-conversion
```

```{toctree}
:hidden:
:caption: Science and implementation

architecture
scientific-guide
validation
fortran-comparison
```

```{toctree}
:hidden:
:caption: Reference

reference/index
```

```{toctree}
:hidden:
:caption: Development

documentation
releasing
```
