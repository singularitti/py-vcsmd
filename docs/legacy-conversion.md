# Legacy conversion

`vcsmd.compat` is the one-way boundary for the historical fixed-format
VCSMD files. The numerical core does not import this adapter. Historical
filenames and two-letter calculation modes are documented here because they
are recognized only while importing old runs.

## Artifacts

`convert_legacy_run(source, destination)` creates a dedicated conversion folder
with this layout:

```text
inputs/
  legacy_source/                 copied raw legacy files
  conversion_manifest.json       relative names, sizes, and SHA-256 hashes
outputs/
  configuration.json             validated native configuration, when inp is usable
  observables.csv                values combined by step_index
  cell_history.csv               source cell moduli and angles, when available
  partial_state.json             structured but incomplete io state, when available
  conversion_report.json         available, missing, malformed, and limited data
README.md                         purpose, provenance, and completion notes
```

Raw files are copied byte-for-byte under `inputs/legacy_source/`; the
manifest records source and converter hashes without machine-specific absolute paths. The
inactive `sip` potential table may be preserved there for archival purposes,
but it is excluded from the computational configuration.

## Imported schemas

The old `inp` reader recognizes `md`, `mm`, `cd`, `cm`, `nd`, `nm`, `sd`, and
`sm` and maps them to the native fixed, cell, metric, and strain modes. It
converts the historical pressure, mass, lattice, and Lennard–Jones parameter
conventions into the native unit records. Fresh inputs return only the strict
native configuration mapping. A `c` continuation imports geometry from
`io` and starts with thermal initialization; it is explicitly not an exact
restart. A historical `f` follow input is refused as a runnable native
configuration, while any readable `io` state remains available as
`partial_state.json`.

The `e` tape supplies potential energy, total kinetic energy (including the
cell contribution), total energy, and `pv`. The `eal` tape supplies atomic and
cell energy components separately. These overlapping quantities remain
separate in `observables.csv`. The `ave` tape contains running potential and
kinetic averages, `p` contains external/internal/average pressure, `tv`
contains cell volume and the running extended temperature, and `avec` contains
cell moduli and angles. `avec` does not contain cell vectors, so no vector
history is reconstructed. Header-only or incomplete `car` files do not become
trajectories because the historical position writer emitted no usable
positions.

Imported `io` records carry explicit units for cell vectors, cell rates,
acceleration history, masses, and fractional particle fields. They also carry
quality metadata for missing absolute origins, reference geometry, Beeman
history, controller accumulators, and duplicated source columns. They are
review data rather than native checkpoints.

## Evidence and limits

Conversion reports retain missing files, malformed records, invalid or
duplicate steps, and known source limitations. Duplicate steps keep the first
record while all raw source bytes remain archived. The supplied full input
document has been parsed and all eight configurations pass native schema
validation. No complete historical output run was available in this
repository, so conversion of real `e`, `eal`, `ave`, `p`, `avec`, `tv`, `car`,
or `io` output tapes remains unvalidated against a complete fixture.
