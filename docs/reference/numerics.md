# Advanced numerical modules

The numerical modules operate on normalized float64 arrays. They are useful
for diagnostics, custom integrations, and understanding the implementation;
ordinary callers should start with the quantity-bearing API.

## Normalized state

```{eval-rst}
.. autoclass:: vcsmd.models.SimulationState
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.RunningAverages
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.Observables
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.SimulationEvent
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.models.StepResult
   :members:
   :undoc-members:
   :show-inheritance:
```

## Geometry

Cell matrices store cell vectors as columns. Fractional particle arrays use one
particle per row, and Cartesian positions are `fractional_positions @ cell.T`.

```{eval-rst}
.. autoclass:: vcsmd.geometry.CellGeometry
   :members:
   :show-inheritance:
```

```{eval-rst}
.. autofunction:: vcsmd.geometry.cell_geometry
```

```{eval-rst}
.. autofunction:: vcsmd.geometry.cell_parameters
```

```{eval-rst}
.. autofunction:: vcsmd.geometry.replicate_cell
```

```{eval-rst}
.. autofunction:: vcsmd.geometry.validate_image_range
```

## Normalized dynamics

These functions operate on `vcsmd.models.NumericalModel` and normalized
states. They are the advanced counterpart of the quantity-bearing functions
in the public API.

```{eval-rst}
.. autofunction:: vcsmd.dynamics.initialize
```

```{eval-rst}
.. autofunction:: vcsmd.dynamics.step
```

```{eval-rst}
.. autofunction:: vcsmd.dynamics.simulate
```

## Lennard–Jones evaluation

```{eval-rst}
.. autoclass:: vcsmd.potentials.ForceResult
   :members:
   :show-inheritance:
```

```{eval-rst}
.. autofunction:: vcsmd.potentials.evaluate_lennard_jones
```

The evaluator uses an unshifted, inclusive cutoff and returns potential energy,
Cartesian forces, and configurational virial. All particles share the global
parameter pair from the normalized model.

## Unit constants

```{eval-rst}
.. autodata:: vcsmd.units.ureg
   :no-value:
```

```{eval-rst}
.. autodata:: vcsmd.units.Quantity
   :no-value:
```

The normalized solver uses bohr, Rydberg, `rydberg_time`, and the corresponding
mass units. Unit conversion at the application boundary is handled by
`vcsmd.config` and `vcsmd.api`.
