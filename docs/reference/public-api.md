# Public simulation API

The user-facing result and evolution types are implemented in `vcsmd.simulation`
and re-exported by `vcsmd`. `vcsmd.SimulationState`,
`vcsmd.Observables`, and `vcsmd.StepResult` are the quantity-bearing classes;
the similarly named classes in `vcsmd.models` are normalized numerical types.

## State and results

```{eval-rst}
.. autoclass:: vcsmd.SimulationState
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.Observables
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.StepResult
   :members:
   :undoc-members:
   :show-inheritance:
```

`StepResult.observables` describes the completed step before temperature
rescaling or minimization quenching. `StepResult.state` is the state after
those actions and is the state to pass to the next step or save.

## Evolution functions

```{eval-rst}
.. autofunction:: vcsmd.initialize
```

```{eval-rst}
.. autofunction:: vcsmd.step
```

```{eval-rst}
.. autofunction:: vcsmd.simulate
```

```{eval-rst}
.. autofunction:: vcsmd.with_units
```

For interactive use, `from vcsmd import *` is supported and follows
`vcsmd.__all__`; explicit imports remain preferable in maintained code.
Configuration loading/saving, checkpoints, and application execution are also
exported from the package root. Their signatures are documented on the
[native I/O](io.md) and [execution](execution.md) pages.
