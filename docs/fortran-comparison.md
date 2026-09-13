# Executed original Fortran comparison

The unmodified `celq.f` was compiled with GNU Fortran 16.2.0 using
`-std=legacy -O2`, without fast-math or source patches. All eight original
Fortran-format inputs completed their full **9,200 integration steps**, with
finite output and contiguous step indices. The reference source, `sip`, input
records, binary, compiler log and all original output tapes are preserved.

The Python comparison reuses the eight complete YAML runs in the
[input-format validation](validation.md), with copied inputs, source, CSV
streams and checkpoints. All 18 Python source modules in each saved run are
byte-identical to the current modules, and each configuration matches the
corresponding original Fortran input after normalization.

## Method and tolerances

`scripts/compare_fortran.py` reads the original fixed-width output independently
of the legacy output converter. Measurements are aligned at `step * input dt`,
after the corrected move and before temperature rescaling or quenching. In
particular, `tv` contains **running extended temperature**, so it is compared
to that Python field. It is not instantaneous atomic temperature. Fortran's
`d12.5` records retain five significant digits; 32 full-stream consistency
checks pass using propagated half-print-unit rounding budgets for energy
decompositions and the pressure-volume product.

The engineering bands were chosen before running Fortran: 5% for energies and
volume, 10% for pressure and temperature, 2% for cell lengths, and 1 degree for
angles. Absolute floors are `N * 1e-8 Ry` for total-cell energy fields,
`1e-8 Ry/bohr^3` for pressure, 1 K for temperature, `1e-5 bohr^3` for volume,
and `1e-5 bohr` for lengths. The allowed difference is the larger of the relative
band and absolute floor. Relative scales use the mean absolute Fortran value
over the comparison window, avoiding unstable percentages near zero pressure.

Dynamic runs compare means over the last half of each full trajectory;
minimizations compare the final step. The last-half window is descriptive,
not a demonstrated equilibrium interval. The table gives absolute percentage
differences between those selected summaries:

| Input | Comparison window | Atomic potential energy | Total energy | Volume |
| --- | --- | ---: | ---: | ---: |
| 1 | Steps 501–1,000 | 0.000613% | 0.000273% | 0.000807% |
| 2 | Steps 1,001–2,000 | 0.002395% | 0.004545% | 0.009387% |
| 3 | Steps 1,001–2,000 | 0.001441% | 0.002951% | 0.008372% |
| 4 | Steps 1,001–2,000 | 0.002412% | 0.004305% | 0.009361% |
| 5 | Final step 100 | 0.001390% | 0.001390% | 0.0000209% |
| 6 | Final step 100 | 0.001583% | 0.001583% | 0.000488% |
| 7 | Steps 501–1,000 | 0.091268% | 0.437145% | 0% |
| 8 | Steps 501–1,000 | 2.552561% | 1.151340% | 0% |

All **92 selected summary checks** are inside their declared bands or absolute
floors. These are correlated checks, including some invariant or zero-valued
quantities, not 92 independent demonstrations of equivalence. Sensitivity at
1%, 5% and 10% is recorded separately in the CSV, retaining the same absolute
floors. Input 5's final running temperature differs by 0.4734 K (10.79%): it
passes the 1 K floor, not the relative 10% band. Final minimized pressures are
close to zero, so their absolute differences are more meaningful than relative
percentages.

## Final structures and limits

An independent postprocessing script, `scripts/compare_fortran_structures.py`,
compares the final Fortran `io` geometry to Python checkpoints for inputs 1–6.
It accounts for cell-vector column ordering and lattice scale, aligns the
fractional origin to the first atom, and compares relative positions modulo
integer cell translations. It skips the damaged legacy cell derivative fields.

The maximum relative cell-matrix difference for dynamic inputs 1–4 is 0.122%.
For the two minimizations it is about `2.8e-8` (0.0000028%), with maximum
cell-vector length differences below `3.3e-7 bohr`. Relative fractional atomic
positions in inputs 1, 2, 5 and 6 differ by at most `4.8e-13`; single-atom
inputs 3 and 4 have no internal relative atomic structure to compare. This
uses saved final states only; no additional trajectory is integrated.

The thermal cases use different random streams: Fortran keeps its hard-coded
legacy seed and Python retains NumPy seed 119. Their individual trajectories
diverge, even though the selected summaries agree. Full-trajectory paired RMSE,
maximum errors, standard deviations and separate third/fourth-quarter means
are retained; no independent-sample confidence interval is claimed. Fortran's
`car` file contains only a header, so particle-by-particle trajectory comparisons
are unavailable.

The original larger timesteps produce measurable energy drift in both codes:
for input 3, the endpoint changes are about -6.30% (Python) and -6.32%
(Fortran). Agreement with this reference is not a claim of energy conservation,
timestep convergence or validation of modes absent from the eight examples.

## Preserved artifacts and reproduction

Artifacts are retained in
`runs/vcsmd-comparison-celq-eight-inputs-results-2026-09-12-21-attempt-01/`:

- `README.md` records creation time, purpose, comparison rules and provenance.
- `inputs/` contains copied reference source and inputs, complete saved Python
  runs, the executed comparison driver, the final analysis scripts, dependency
  versions and SHA-256 manifests. The original execution protocol is retained;
  later clarifications correct its output-timing/precision description without
  changing the bands, windows or simulation outputs.
- `outputs/build/` contains `celq` and its compiler log.
- `outputs/fortran/` preserves every full original output tape, execution log
  and normalized CSV. Each run documents its `inp`/`sip` symlinks to copied inputs.
- `outputs/metrics.csv`, `run_summary.csv`, `consistency_checks.csv`,
  `energy_conservation.csv` and `final_structures.csv` contain quantitative results.
- `outputs/REPORT.md` and the eight-page `diagnostics.pdf` present the results;
  `outputs/output_hashes.csv` records their hashes.

Reproduce all eight full reference trajectories and their comparison after
installing the project's dev dependencies:

```bash
uv run python scripts/compare_fortran.py \
  --legacy-directory ../vcsmd \
  --python-suite runs/vcsmd-validation-yaml-legacy-inputs-results-2026-09-12-21-attempt-01
```

Regenerate analysis and the final-structure check from the preserved outputs,
without running trajectories:

```bash
uv run python scripts/compare_fortran_structures.py --comparison-suite \
  runs/vcsmd-comparison-celq-eight-inputs-results-2026-09-12-21-attempt-01
uv run python scripts/compare_fortran.py --analyze \
  runs/vcsmd-comparison-celq-eight-inputs-results-2026-09-12-21-attempt-01
```

Running the main analysis last refreshes the manifests for all generated files.
