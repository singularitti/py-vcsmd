# Implementation validation

Validation followed feature implementation. No reduced, synthetic, or altered simulation workloads were created. The original eight examples were converted with their original particle counts, physical parameters, and total step counts. The default local NumPy seed was 119; reproducing the historical random stream was deliberately outside scope.

## Full original workloads

The successful suite is preserved under `runs/vcsmd-validation-eight-inputs-results-2026-09-12-05-attempt-03/`. Every workload completed with finite observables:

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

The first suite attempt, in the adjacent `attempt-01` directory, is retained. Three workloads completed and five encountered a checkpoint validation error when a negative fractional coordinate very close to zero wrapped to exactly 1 in floating-point arithmetic. The fix maps such remainders to zero. The second suite completed all eight workloads. A final serialization check then found that converting already-normalized quantities through SI could move a mass by one rounding unit on each save/reload. The fix converts quantities directly within the native registry. The final third suite repeated every original workload with that fix and the bounded-memory periodic evaluator. All three attempts retain the same requested physical parameters; earlier outputs and snapshots remain preserved.

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

The final periodic-batch implementation was compared with the preserved full-run implementation at all eight original initial structures and all eight saved final states. That comparison is preserved in the second attempt and checks exact energy, force, and virial equality before the final full-suite repetition. Its CSV is copied into the third attempt with provenance; the original comparison script and its baseline source remain in the second attempt.

## Packaging and adapters

Installation, all module imports, CLI help, equivalent configuration normalization, and the public quantity-bearing checkpoint roundtrip were checked on Python 3.10 and 3.14 with compatible dependency resolution. Exact interpreter and dependency versions are recorded in `outputs/package_checks_python310.csv` and `outputs/package_checks_python314.csv`. The public checkpoint check used the complete final state from input 2. CLI input conversion used an unchanged original input block. Temporary serialization files were removed.

The available historical source directory contains no complete historical output tapes. Conversion correctly reported missing observables, trajectory, and final restart data, and all eight original input blocks parsed through native schema validation. Historical output parsing received structural checks, but conversion of a real complete historical run remains unverified. Imported values are not claimed to be recomputed or scientifically corrected.

## Preserved artifacts and reproduction

Within the successful suite folder:

- `inputs/` contains the copied original input collection and its hash, converted full configurations, configuration-format variants, and exact validation scripts. Additional checks are identified separately from the original suite snapshot.
- `outputs/summary.csv`, `checks.csv`, and `numerical_checks.csv` record workload outcomes and checks. `force_batch_checks.csv` records the final batch comparison. Package checks have separate CSV files per interpreter.
- `outputs/diagnostics.pdf` contains the full-run plots.
- `outputs/runs/` contains native run directories with CSV observables, cell histories, particle trajectories, controller events, periodic and final NPZ checkpoints, source copies, dependency metadata, and run status. The complete continuation output is retained.
- `outputs/validation_metadata.json` records coverage and suite outcome. Conversion directories contain raw copied historical artifacts and missing-data reports.

The completed simulation artifacts and their original source snapshots are retained. Subsequent documentation and converter-report fixes are recorded as post-validation changes; original snapshots are not overwritten. Source hashes, rather than a Git revision alone, identify the uncommitted implementation used for each run. Run artifacts are ignored by Git and remain in the workspace.

To repeat the full acceptance workflow after installation:

```bash
uv run python scripts/validate_full_examples.py --legacy-directory ../vcsmd
uv run python scripts/check_numerics.py runs/<full-validation-directory>
```

Five implemented modes have no supplied full-trajectory cases: `fixed-min`, `cell-dyn`, `cell-min`, `strain-dyn`, and `strain-min`. Their implementation has equation/source review but must not be described as fully trajectory-validated. No Fortran compiler/reference execution was available, so direct numerical parity is unverified. The attached paper's electronic-structure calculations and tensor external stress remain outside the scientific scope.
