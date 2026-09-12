# Configuration and normalized models

`vcsmd.config` contains the unit-aware domain objects used to construct a
simulation. `prepare` converts them into immutable normalized objects for the
solver. Native JSON, YAML, and TOML files are parsed into these same objects by
the functions in the [I/O reference](io.md).

## Unit-aware configuration

```{eval-rst}
.. autoclass:: vcsmd.config.Structure
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.config.LennardJones
   :members:
   :undoc-members:
   :show-inheritance:
```

All particles use one global Lennard–Jones parameter pair; `species` labels do
not select separate interactions or mixing rules.

```{eval-rst}
.. autoclass:: vcsmd.config.SimulationConfig
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autofunction:: vcsmd.config.prepare
```

## Enums and controller settings

These definitions are also re-exported from `vcsmd`:

```{eval-rst}
.. autoclass:: vcsmd.models.SimulationMode
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.InitializationMode
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.EventKind
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.TemperatureControl
   :members:
   :undoc-members:
   :show-inheritance:
```

Variable-cell modes require `cell_inertia`. Ordinary cell modes use mass units;
modified-metric and reference-strain modes use mass/length⁴. Provided particle
velocities must satisfy the zero mass-weighted center-of-mass constraint.

## Normalized model inputs

These immutable types are the advanced numerical boundary returned by
`prepare`. Their arrays use internal units and are read-only.

```{eval-rst}
.. autoclass:: vcsmd.models.NumericalModel
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.InitialConditions
   :members:
   :undoc-members:
   :show-inheritance:
```
