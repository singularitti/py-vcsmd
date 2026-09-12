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
