# Implementation validation

Validation followed feature implementation. No reduced, synthetic, or altered simulation workloads were created. The original eight examples were converted with their original particle counts, physical parameters, and total step counts. The default local NumPy seed was 119; reproducing the historical random stream was deliberately outside scope.

## Shipped YAML and individual legacy inputs

The earlier checks below established serialization equivalence and parsed all
eight legacy blocks, but did not execute every shipped YAML file separately.
The examples now include all eight YAML configurations and eight individual
Fortran-format `.inp` files. The latter preserve the original input records,
removing leading blank collection separators so the title is the first record.

The dedicated format validation completed **16 full runs and 89 checks**:

| Input path | Full workloads completed | Integration steps |
| --- | ---: | ---: |
| Shipped YAML loaded through `load_config` | 8/8 | 9,200 |
| Individual legacy inputs converted by the `import-legacy` CLI, then loaded from YAML | 8/8 | 9,200 |

All original particle counts, physical parameters, and requested step counts
were preserved. Historical inputs contain no random seed; the importer uses
the native default seed 119, matching the shipped YAML configurations.
Every run completed with finite observables at every step. For each
workload, the YAML and imported legacy configurations matched the canonical
configuration parsed from the source collection and the shipped TOML file.
The first JSON file also matched. Each YAML/legacy pair produced exactly equal
final checkpoint arrays and byte-identical complete CSV output streams.
These are internal consistency checks of the Python import-and-run workflows.
The subsequent [original Fortran comparison](fortran-comparison.md) compiled
unmodified `celq.f` and completed all eight full reference workloads. It compares
energies, pressure, temperature, volume and final structures using declared
tolerances, with the full results and limitations documented separately.

The local suite is retained under
`runs/vcsmd-validation-yaml-legacy-inputs-results-2026-09-12-21-attempt-01/`:

- `README.md` records creation time, purpose, input provenance, contents, and outcome.
- `inputs/examples/` contains copied source examples, and `inputs/source_hashes.csv`
  records their hashes and the copied validation script's hash.
- `inputs/imported/` contains the native YAML files produced by the CLI.
- `outputs/summary.csv` records all 16 full runs; `outputs/checks.csv` records
  all 89 checks, and `outputs/import-*.log` records CLI import results.
- `outputs/runs/` retains the complete CSV results, periodic and final checkpoints,
  and per-run input/source snapshots with dependency versions and hashes.

Run artifacts are ignored by Git. Reproduce this check from the repository
root after installation:

```bash
uv run python scripts/validate_input_formats.py
```

The script executes all 18,400 steps; it does not create reduced test cases.
The earlier implementation-validation results below describe the original
suite. Its historical run folder is not present in this checkout; the new
format-validation folder above contains the evidence from this verification.

## Original Fortran reference comparison

The [executed original Fortran comparison](fortran-comparison.md) compiled
unmodified `celq.f` and completed all eight full reference workloads. It
compares energies, pressure, temperature, volume and final structures using
declared tolerances, with preserved inputs, outputs and reproducible scripts.

## Full original workloads

The original implementation report recorded the successful suite under `runs/vcsmd-validation-eight-inputs-results-2026-09-12-05-attempt-03/`. Its reported results follow; every workload completed with finite observables:

| Original input | Mode | Particles | Requested and completed steps |
| --- | --- | ---: | ---: |
| 1 | `metric-dyn` | 4 | 1,000 |
| 2 | `metric-dyn` | 4 | 2,000 |
| 3 | `metric-dyn` | 1 | 2,000 |
| 4 | `metric-dyn` | 1 | 2,000 |
| 5 | `metric-min` | 4 | 100 |
| 6 | `metric-min` | 2 | 100 |
| 7 | `fixed-dyn` | 32 | 1,000 |
| 8 | `fixed-dyn` | 32 | 1,000 |

This totals 9,200 steps. Input 2 was also continued from its step-1,000 checkpoint through its original final step 2,000. Every field of the complete final numerical state was exactly equal to the uninterrupted run, including accelerations, reference geometry, controller counter, and running sums. This verifies continuation within that environment; it does not guarantee bitwise identity across different numerical libraries or hardware.

The first suite attempt completed three workloads; its failed and intermediate artifacts were removed during the requested cleanup. Five workloads encountered a checkpoint validation error when a negative fractional coordinate very close to zero wrapped to exactly 1 in floating-point arithmetic. The fix maps such remainders to zero. The second suite completed all eight workloads. A final serialization check then found that converting already-normalized quantities through SI could move a mass by one rounding unit on each save/reload. The fix converts quantities directly within the native registry. The final third suite repeated every original workload with that fix and the bounded-memory periodic evaluator. All three attempts used the same requested physical parameters. Only the final third suite is retained; failed and superseded artifacts were removed at the user's request.

Completion is not a convergence assertion. The relative total-energy ranges for inputs 1–4 were approximately 4.152%, 0.643%, 6.953%, and 0.643%, respectively. These original timesteps were preserved; no timestep study was performed. Inputs 5–6 quench velocities and inputs 7–8 apply temperature rescaling, so their whole-run energy ranges are not conservation measures. Energy, pressure, volume, and both atomic and running extended temperatures are plotted in the preserved PDF.

## Post-feature checks

All 18 full-suite checks and all 17 numerical checks passed. Numerical derivative probes evaluated complete supplied systems from their saved states; no extra trajectory was integrated.

| Identity or property | Observed maximum absolute discrepancy |
| --- | ---: |
| Full 32-particle force versus negative energy gradient | 2.34 × 10⁻¹² in normalized units |
| Periodic virial versus negative cell-deformation energy derivative | 8.68 × 10⁻¹³ in normalized units |
| Modified-metric cell acceleration versus differentiated Lagrangian | 1.61 × 10⁻¹⁸ in normalized units |
| Global momentum after thermal initialization | 2.58 × 10⁻¹⁴ in normalized units |
| Initialized temperature versus requested temperature | 0 K |

Other checks covered deterministic initialization, unchanged supplied arrays, isolation from NumPy's global random state, read-only published arrays, correct cofactor orientation, periodic force balance, formulation-specific inertia dimensions, half-open wrapping, invalid average counters, and absence of I/O, compatibility, and Pint imports in the normalized numerical modules. JSON, safe YAML, and TOML produced equivalent native models.

The final periodic-batch implementation was compared with the preserved full-run implementation at all eight original initial structures and all eight saved final states. That comparison established exact energy, force, and virial equality before the final full-suite repetition. Its CSV remains as a historical result in the retained final suite. The superseded implementation and comparison script were removed during cleanup, so that particular before/after comparison is no longer reproducible from the retained artifacts alone. The final full trajectories and their exact executed implementation are retained.

## Packaging and adapters

Installation, all module imports, CLI help, equivalent configuration normalization, and the public quantity-bearing checkpoint roundtrip were checked on Python 3.10 and 3.14 with compatible dependency resolution. Exact interpreter and dependency versions are recorded in `outputs/package_checks_python310.csv` and `outputs/package_checks_python314.csv`. The public checkpoint check used the complete final state from input 2. CLI input conversion used an unchanged original input block. Temporary serialization files were removed.

The available historical source directory contains no complete historical output tapes. Conversion correctly reported missing observables, trajectory, and final restart data, and all eight original input blocks parsed through native schema validation. Historical output parsing received structural checks, but conversion of a real complete historical run remains unverified. Imported values are not claimed to be recomputed or scientifically corrected.

## Preserved artifacts and reproduction

Within the successful suite folder:

- `inputs/` contains the copied original input collection and its hash, converted full configurations, configuration-format variants, and exact validation scripts. The nine byte-identical executed-source snapshots were consolidated into one shared `inputs/source/` copy, with dependency declarations and original source hashes. Each child run records its relative location.
- `outputs/summary.csv`, `checks.csv`, and `numerical_checks.csv` record workload outcomes and checks. `force_batch_checks.csv` records the final batch comparison. Package checks have separate CSV files per interpreter.
- `outputs/diagnostics.pdf` contains the full-run plots.
- `outputs/runs/` contains native run directories with CSV observables, cell histories, particle trajectories, controller events, periodic and final NPZ checkpoints, input configurations, and run status. The complete continuation output is retained.
- `outputs/validation_metadata.json` records coverage and suite outcome. Conversion directories contain raw copied historical artifacts and missing-data reports.

The final completed simulation artifacts and one exact executed-source snapshot are retained. Cleanup removed the failed first attempt, superseded second attempt, old conversion report, redundant source copies, and disposable build archives. `outputs/cleanup_inventory.csv` lists the removals, while `inputs/cleanup_provenance.json` records the source consolidation. Later source docstrings and the architecture walkthrough do not change numerical behavior. Source hashes identify the implementation used for each run. Run artifacts are ignored by Git and remain in the workspace.

To repeat the full acceptance workflow after installation:

```bash
uv run python scripts/validate_full_examples.py --legacy-directory ../vcsmd
uv run python scripts/check_numerics.py runs/<full-validation-directory>
```

Five implemented modes have no supplied full-trajectory cases: `fixed-min`, `cell-dyn`, `cell-min`, `strain-dyn`, and `strain-min`. Their implementation has equation/source review but must not be described as fully trajectory-validated. The [executed Fortran comparison](fortran-comparison.md) covers the eight supplied workloads in the other three modes, using declared tolerances rather than exact equality. The attached paper's electronic-structure calculations and tensor external stress remain outside the scientific scope.
