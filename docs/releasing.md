# Publishing a release

The build workflow runs on pushes and pull requests targeting `main`. It builds
the source distribution and wheel on Python 3.9 through 3.14, checks distribution
metadata, installs the wheel and its dependencies, imports the installed modules,
and checks both CLI entry points. These checks do not run simulations.

## One-time PyPI setup

Configure a GitHub Actions Trusted Publisher on PyPI for this repository and the
workflow filename `python-publish.yml`. For a new PyPI project, configure a pending
publisher for the project name declared in `pyproject.toml`.

The workflow uses OpenID Connect with `id-token: write`; no PyPI API token is
required. It currently does not use a GitHub environment, so leave the environment
field empty in the publisher configuration. If you configure an environment on
PyPI, add the matching `environment` to the workflow's `publish` job as well.

## Release steps

1. Choose an unused package version and update both `pyproject.toml` and
   `src/vcsmd/__init__.py`. Run `uv lock` to refresh the lockfile. The initial
   `0.1.0` version can be kept if it has not already been published.
2. Commit and push the changes, including the workflows, and confirm that the
   build matrix passes on `main`.
3. Create a tag at that commit and publish a GitHub release for the tag. The
   publishing workflow checks out the release tag, builds and validates the
   distributions, and uploads them to PyPI.

For a manual publication, run **Upload Python Package** from the Actions tab and
provide the tag or commit in the `ref` input. The workflow must be present on the
default branch to enable manual runs. Manual runs also upload to production PyPI.

A Git tag does not change the package version. PyPI does not allow published
distribution filenames to be reused, so publishing another tag or rerunning a
successful publication requires a new version for new distribution files.
