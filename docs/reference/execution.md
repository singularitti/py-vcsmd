# Execution and run folders

`vcsmd.execution` is the file-writing application boundary. The numerical API
does not access files or emit progress output. `run` creates a documented run
folder with configuration provenance, CSV streams, checkpoints, and metadata;
`resume` continues a native checkpoint for additional steps.

```{eval-rst}
.. autoclass:: vcsmd.execution.RunStatus
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autoclass:: vcsmd.execution.RunReport
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. autofunction:: vcsmd.execution.run
```

```{eval-rst}
.. autofunction:: vcsmd.execution.resume
```

`resume` preserves acceleration history, reference geometry, running
accumulators, and controller state from the checkpoint. Exact continuation is
expected within the same numerical software environment; changes to numerical
libraries, hardware, or code can change floating-point rounding.

The installed `vcsmd` command exposes `run`, `resume`, `import-legacy`, and
`convert-legacy`. See the guides for command examples and the generated output
layout.
