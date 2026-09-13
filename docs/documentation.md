# Building the documentation

The documentation uses [Sphinx](https://www.sphinx-doc.org/),
[MyST Markdown](https://myst-parser.readthedocs.io/), and the
[Furo theme](https://furo.readthedocs.io/). Sphinx reads the Python docstrings to
generate API references. MyST renders the Markdown guides, mathematics, and
Mermaid architecture diagram.

## Install the tools

Run these commands from the repository root with **Python 3.12 or newer**:

```bash
uv sync --group docs --no-default-groups
```

If your selected Python is older, select a compatible interpreter explicitly:

```bash
uv sync --group docs --no-default-groups --python 3.14
```

The `docs` dependency group is separate from runtime dependencies. The package
declares Python 3.10+ support. `uv sync` creates a local `uv.lock` automatically,
but the lockfile is ignored by Git and is not required in a checkout. CI
resolves the dependency ranges in `pyproject.toml` on each fresh run. Builds
can therefore pick up newer compatible releases over time.

## Build HTML

```bash
uv run --group docs --no-default-groups sphinx-build -W --keep-going -b html docs docs/_build/html
```

The site is written to `docs/_build/html/index.html`. `-W` makes warnings fail the
build, including broken document references, invalid directives, and duplicate
API targets. `--keep-going` reports all encountered warnings in the same build.
The generated output is ignored by Git.

The build imports package modules for API introspection. It does not execute
the displayed simulation code, run numerical checks, or create simulation run
folders. Example downloads are copied from `examples/`; simulation results in
`runs/` are not bundled into the website.

## Preview while editing

```bash
uv run --group docs --no-default-groups sphinx-autobuild --watch src --port 8000 docs docs/_build/html
```

Open `http://localhost:8000` in a browser. Changes to documentation and package
source trigger a rebuild. Stop the preview with **Ctrl+C**.

The theme and search are bundled with the site. MathJax and Mermaid are loaded
from version-pinned CDN URLs in `docs/conf.py`, so equation and diagram rendering
requires network access in the browser.

## Update the content

- Edit guides in `docs/` and add pages to a toctree in `docs/index.md` or
  `docs/reference/index.md`.
- Edit Python docstrings to update API descriptions. Reference pages select
  public objects with Sphinx `autoclass` and `autofunction` directives inside
  MyST `eval-rst` fences. Keep these fences: they ensure Sphinx 9's generated
  reStructuredText is parsed into API objects instead of displayed as text.
- Use `$...$` for inline mathematics and `$$...$$` for display mathematics.
- Use fenced `mermaid` blocks for diagrams. The existing architecture guide
  remains readable as Markdown in the repository.
- Keep example listings tied to their original files with `literalinclude`
  and download roles. Preserve full workloads and their provenance.
- Run the warning-free HTML build before submitting documentation changes.

The **Build Documentation** workflow builds the same site on pushes and pull
requests, checks that API objects were rendered, and uploads the generated
HTML as a downloadable Actions artifact. It can also be run manually.
Publication to a hosting service is a separate step.
