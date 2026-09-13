# API reference

This reference is organized around the supported application boundary. The
top-level `vcsmd` package re-exports the most common names from the modules
below.

```{toctree}
:maxdepth: 2

public-api
configuration
execution
io
numerics
compatibility
```

The public API uses Pint quantities and read-only arrays. The normalized model
and state types in `vcsmd.models` use float64 values in the package's internal
units and are intended for the numerical kernels, persistence, and advanced
callers.

The supported application surface is re-exported from the top-level `vcsmd`
package. Its `__all__` defines the names imported by `from vcsmd import *`, which
is convenient for interactive exploration. Maintained application and library
code should generally use explicit imports so dependencies remain visible.
