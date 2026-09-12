# Public simulation API

The user-facing result and evolution types live in `vcsmd.api` and are
re-exported by `vcsmd` for convenience. `vcsmd.SimulationState`,
`vcsmd.Observables`, and `vcsmd.StepResult` are the quantity-bearing classes;
the similarly named classes in `vcsmd.models` are normalized numerical types.

## State and results

```{eval-rst}
.. autoclass:: vcsmd.api.SimulationState
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.api.Observables
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.api.StepResult
   :members:
   :undoc-members:
   :show-inheritance:
```

`StepResult.observables` describes the completed step before temperature
rescaling or minimization quenching. `StepResult.state` is the state after
those actions and is the state to pass to the next step or save.

## Evolution functions

```{eval-rst}
.. autofunction:: vcsmd.api.initialize
```

```{eval-rst}
.. autofunction:: vcsmd.api.step
```

```{eval-rst}
.. autofunction:: vcsmd.api.simulate
```

```{eval-rst}
.. autofunction:: vcsmd.api.with_units
```

The equivalent top-level imports are `vcsmd.initialize`, `vcsmd.step`,
`vcsmd.simulate`, and `vcsmd.with_units`.
