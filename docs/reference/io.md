# Native configuration and checkpoint I/O

The `vcsmd.io` namespace is the recommended facade over the format-specific
modules.

## Configuration mappings and files

```{eval-rst}
.. autofunction:: vcsmd.io.config.config_from_mapping
```

```{eval-rst}
.. autofunction:: vcsmd.io.config.config_to_mapping
```

```{eval-rst}
.. autofunction:: vcsmd.io.config.load_config
```

```{eval-rst}
.. autofunction:: vcsmd.io.config.save_config
```

`load_config` selects JSON, YAML, or TOML from the path suffix. All three
formats use the same schema and validation pipeline. Dimensional scalars are
records with `value` and `unit`; dimensional arrays use `values` and `unit`.

## Checkpoints

```{eval-rst}
.. autofunction:: vcsmd.io.load_checkpoint
```

```{eval-rst}
.. autofunction:: vcsmd.io.save_checkpoint
```

The facade wraps a loaded numerical state as `vcsmd.api.SimulationState` and
requires that public state when saving. Native archives use named numeric
arrays and do not use Python pickle objects.

The lower-level checkpoint implementation is available for advanced numerical
code:

```{eval-rst}
.. autofunction:: vcsmd.io.checkpoint.load_checkpoint
```

```{eval-rst}
.. autofunction:: vcsmd.io.checkpoint.save_checkpoint
```
