# Native configuration and checkpoint I/O

The top-level `vcsmd` namespace is the recommended facade over the
format-specific I/O modules. The implementation remains organized under
`vcsmd.io`.

## Configuration mappings and files

```{eval-rst}
.. autofunction:: vcsmd.config_from_mapping
```

```{eval-rst}
.. autofunction:: vcsmd.config_to_mapping
```

```{eval-rst}
.. autofunction:: vcsmd.load_config
```

```{eval-rst}
.. autofunction:: vcsmd.save_config
```

`load_config` selects JSON, YAML, or TOML from the path suffix. All three
formats use the same schema and validation pipeline. Dimensional scalars are
records with `value` and `unit`; dimensional arrays use `values` and `unit`.

## Checkpoints

```{eval-rst}
.. autofunction:: vcsmd.load_checkpoint
```

```{eval-rst}
.. autofunction:: vcsmd.save_checkpoint
```

The facade wraps a loaded numerical state as `vcsmd.SimulationState` and
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
