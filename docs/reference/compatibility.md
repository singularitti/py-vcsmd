# Legacy compatibility

`vcsmd.compat` is an isolated, one-way adapter for historical fixed-format
inputs and output directories. Imported data preserves available source fields
and provenance. Missing or damaged restart history remains partial and is not
represented as an exact native checkpoint.

```{eval-rst}
.. autoclass:: vcsmd.compat.LegacyFormatError
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.compat.ConversionReport
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autofunction:: vcsmd.compat.parse_legacy_input
```

```{eval-rst}
.. autofunction:: vcsmd.compat.import_legacy_input
```

```{eval-rst}
.. autofunction:: vcsmd.compat.split_legacy_examples
```

```{eval-rst}
.. autofunction:: vcsmd.compat.convert_legacy_run
```

The compatibility guide describes recognized historical files, imported CSV
schemas, and conversion limitations.
