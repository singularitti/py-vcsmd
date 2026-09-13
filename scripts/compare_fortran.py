"""Compile original celq.f and compare all eight full workloads to saved Python runs.

No trajectory is shortened and no simulation parameters or random seeds change.
Use --analyze RUN_DIRECTORY to regenerate analysis without integrating again.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from vcsmd import config_from_mapping, config_to_mapping, load_config
from vcsmd.compat import parse_legacy_input

PROTOCOL = """# Full original Fortran versus Python comparison

The reference is the copied, unmodified celq.f, compiled with gfortran in legacy
mode at -O2 without fast-math. All eight original workloads retain their input
records, particle counts, physical parameters and full step counts. Fortran
uses its hard-coded legacy seed (-119); Python retains its configured seed
(119). The thermal random streams differ; neither initialization is changed.
The Python side reuses complete saved YAML runs; their input configuration must
match the corresponding Fortran input after normalization. Every saved Python
run, including its executed source, input, checkpoints and raw outputs, is copied.

Before running the reference, the practical agreement bands are set to 5% for
energies and volume, 10% for pressure and temperature, 2% for cell lengths, and
1 degree for cell angles. Absolute floors are N * 1e-8 Ry for each total-cell
energy field (including cell kinetic energy), where N is the particle count,
1e-8 Ry/bohr^3 for pressure, 1 K for temperature, 1e-5 bohr^3 for volume,
and 1e-5 bohr for lengths. Relative scales are mean(abs(reference)) in the
comparison window, so pressure zero crossings do not create singular ratios.
The allowed difference is max(relative band * scale, absolute floor).
These bands are engineering screens, not fitted error bars or a proof of parity.

Dynamic runs use last-half means as the primary comparison; full-run means,
first/final values, paired RMSE and maximum absolute errors are also reported.
Minimizations use the final step, without asserting that 100 steps converged.
The last half is a fixed descriptive window, not a demonstrated equilibrium
window. Its two quarter means and standard deviation reveal drift/fluctuations.
Frames are correlated, initialization differs for thermal runs, and there is
only one trajectory per implementation: no confidence interval is claimed.
Sensitivity screens at 1%, 5%, and 10% accompany every non-angle observable.
Those use the same absolute floor and are not changes to the primary bands.

Fortran energy tapes are in Ry (total cell energy), pressure in Ry/bohr^3,
volume in bohr^3, lengths in bohr and angles in degrees. tv temperature is a
running extended temperature, matched to Python's running_extended_temperature.
Fortran calls move (including correction and measurement) before updating its
averages and saving observables; both streams therefore use step * input dt
and describe the completed step before temperature rescaling or quenching.
Fortran's printed d12.5 values (for example, 0.20539D-01) retain five
significant decimal digits. Tape identities use the sum of half-print-unit
rounding uncertainties, with product uncertainties propagated for pressure*volume.

Inputs are copied snapshots under inputs/: original source, sip, legacy inputs,
comparison script and saved Python runs. Symlinks named inp and sip inside each
Fortran output directory refer to those copied inputs, not live originals.
Outputs include the compiler log and binary, all original Fortran tapes and
logs, normalized CSV streams, run_summary.csv, metrics.csv, consistency_checks.csv,
energy_conservation.csv, diagnostics.pdf and REPORT.md. The optional saved-state
postprocessor adds final_structures.csv and copies compare_fortran_structures.py
under inputs; it does not integrate trajectories. source_hashes.csv and
output_hashes.csv record SHA-256 hashes.
The original source and previous Python runs remain in their original locations.
No reduced, temporary or synthetic simulation cases are created.
"""

# metric: Python column, Fortran tape/column, units, relative band, absolute floor
METRICS = {
    "atomic_potential_energy": (
        "atomic_potential_energy [Ry]",
        "eal",
        0,
        "Ry",
        0.05,
        1e-8,
    ),
    "atomic_kinetic_energy": ("atomic_kinetic_energy [Ry]", "eal", 1, "Ry", 0.05, 1e-8),
    "cell_kinetic_energy": ("cell_kinetic_energy [Ry]", "eal", 4, "Ry", 0.05, 1e-8),
    "total_energy": ("total_energy [Ry]", "e", 2, "Ry", 0.05, 1e-8),
    "pressure": ("pressure [Ry/bohr^3]", "p", 1, "Ry/bohr^3", 0.10, 1e-8),
    "volume": ("volume [bohr^3]", "tv", 0, "bohr^3", 0.05, 1e-5),
    "running_extended_temperature": (
        "running_extended_temperature [K]",
        "tv",
        1,
        "K",
        0.10,
        1.0,
    ),
    "a": ("a [bohr]", "avec", 0, "bohr", 0.02, 1e-5),
    "b": ("b [bohr]", "avec", 1, "bohr", 0.02, 1e-5),
    "c": ("c [bohr]", "avec", 2, "bohr", 0.02, 1e-5),
    "gamma": ("gamma [degree]", "avec", 3, "degree", 0.0, 1.0),
    "alpha": ("alpha [degree]", "avec", 4, "degree", 0.0, 1.0),
    "beta": ("beta [degree]", "avec", 5, "degree", 0.0, 1.0),
}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def hashes(root: Path, destination: Path) -> None:
    write_csv(
        destination,
        [
            {
                "file": p.relative_to(root).as_posix(),
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
            for p in sorted(root.rglob("*"))
            if p.is_file() and not p.is_symlink() and p != destination
        ],
    )


def tape(path: Path, columns: int, steps: int) -> np.ndarray:
    """Independently read the original 1x,Nd12.5,i6 writer, including finite checks."""
    lines = path.read_text().splitlines()
    if int(lines[0]) != steps or len(lines) != steps + 1:
        raise ValueError(f"{path.name}: incomplete trajectory")
    result = []
    for number, line in enumerate(lines[1:], 1):
        values = [
            float(line[1 + 12 * i : 1 + 12 * (i + 1)].replace("D", "E"))
            for i in range(columns)
        ]
        if int(line[1 + 12 * columns :]) != number or not np.isfinite(values).all():
            raise ValueError(f"{path.name}: invalid record {number}")
        result.append(values)
    return np.asarray(result)


def print_uncertainty(values: np.ndarray) -> np.ndarray:
    """Half of the last retained digit in the original 0.dddddD+ee format."""
    magnitude = np.maximum(np.abs(values), np.finfo(float).tiny)
    return np.where(values == 0, 0.0, 0.5 * 10.0 ** (np.floor(np.log10(magnitude)) - 4))


def analyze(directory: Path) -> None:
    output = directory / "outputs"
    original_protocol = directory / "inputs/execution_protocol.md"
    if not original_protocol.exists():
        shutil.copyfile(directory / "README.md", original_protocol)
    metadata = original_protocol.read_text().split("\nCreated:", 1)[1]
    (directory / "README.md").write_text(
        PROTOCOL
        + "\nCreated:"
        + metadata
        + "\nAnalysis uses the corrected output-timing description above; the original execution driver and protocol are retained under inputs. Numerical bands and comparison windows are unchanged.\n"
    )
    # Preserve the actual execution driver and the later analysis revision separately.
    analysis_snapshot = directory / "inputs/analyze_fortran.py"
    if Path(__file__).resolve() != analysis_snapshot.resolve():
        shutil.copyfile(__file__, analysis_snapshot)
    write_csv(
        directory / "inputs/analysis_environment.csv",
        [
            {"component": "Python", "version": sys.version.split()[0]},
            {"component": "NumPy", "version": np.__version__},
            {"component": "Matplotlib", "version": matplotlib.__version__},
        ],
    )
    hashes(directory / "inputs", directory / "inputs/source_hashes.csv")
    metrics = []
    conservation = []
    checks = []
    completed = []
    with PdfPages(output / "diagnostics.pdf") as pdf:
        for row in read_csv(output / "run_summary.csv"):
            number, steps, particles = (
                int(row[k]) for k in ("input", "steps", "particles")
            )
            if row["status"] != "complete":
                continue
            reference = output / row["fortran_directory"]
            saved = directory / "inputs/python_runs" / f"input-{number:02d}" / "outputs"
            py = read_csv(saved / "observables.csv")
            cells = read_csv(saved / "cell_history.csv")
            if [int(r["step_index"]) for r in py] != list(range(1, steps + 1)):
                raise ValueError(f"Python input {number}: incomplete trajectory")
            tables = {
                name: tape(reference / name, count, steps)
                for name, count in {"e": 4, "eal": 6, "p": 3, "tv": 2}.items()
            }
            if row["mode"] != "fixed-dyn":
                tables["avec"] = tape(reference / "avec", 6, steps)
            # Budget each rounded field independently, including cancellation.
            e, eal, p, tv = (tables[name] for name in ("e", "eal", "p", "tv"))
            de, deal, dp, dtv = (print_uncertainty(t) for t in (e, eal, p, tv))
            identities = {
                "total_energy_equals_potential_plus_kinetic": (
                    e[:, 2],
                    e[:, 0] + e[:, 1],
                    de[:, 2] + de[:, 0] + de[:, 1],
                ),
                "atomic_plus_cell_kinetic": (
                    e[:, 1],
                    eal[:, 1] + eal[:, 4],
                    de[:, 1] + deal[:, 1] + deal[:, 4],
                ),
                "atomic_plus_external_potential": (
                    e[:, 0],
                    eal[:, 0] + eal[:, 3],
                    de[:, 0] + deal[:, 0] + deal[:, 3],
                ),
                "pressure_volume_product": (
                    e[:, 3],
                    p[:, 1] * tv[:, 0],
                    de[:, 3]
                    + np.abs(p[:, 1]) * dtv[:, 0]
                    + np.abs(tv[:, 0]) * dp[:, 1]
                    + dp[:, 1] * dtv[:, 0],
                ),
            }
            for name, (left, right, uncertainty) in identities.items():
                passed = bool(np.all(np.abs(left - right) <= uncertainty + 1e-14))
                checks.append(
                    {
                        "input": number,
                        "check": name,
                        "passed": passed,
                        "max_absolute_residual": float(np.max(np.abs(left - right))),
                    }
                )
                if not passed:
                    raise ValueError(
                        f"Input {number}: Fortran tape identity {name} failed"
                    )
            time = np.asarray([float(r["time [rydberg_time]"]) for r in py])
            config = load_config(
                directory
                / "inputs/python_runs"
                / f"input-{number:02d}"
                / "inputs/configuration.json"
            )
            expected_time = (
                np.arange(1, steps + 1) * config.timestep.to("rydberg_time").magnitude
            )
            if not np.allclose(time, expected_time, rtol=1e-12, atol=0):
                raise ValueError(
                    f"Input {number}: Python time axis differs from step * input dt"
                )
            if [int(r["step_index"]) for r in cells] != list(range(1, steps + 1)):
                raise ValueError(f"Input {number}: incomplete Python cell history")
            normalized = [
                {"step_index": i + 1, "time [rydberg_time]": time[i]}
                for i in range(steps)
            ]
            plot_data = {}
            for name, (
                column,
                filename,
                index,
                unit,
                relative,
                floor,
            ) in METRICS.items():
                if filename not in tables:
                    continue
                p = np.asarray(
                    [float(r[column]) for r in (cells if filename == "avec" else py)]
                )
                f = tables[filename][:, index]
                if not np.isfinite(p).all():
                    raise ValueError(f"Python input {number}: nonfinite {name}")
                if unit == "Ry":
                    floor *= particles
                for i, value in enumerate(f):
                    normalized[i][f"{name} [{unit}]"] = value
                plot_data[name] = (p, f)
                window = (
                    slice(-1, None)
                    if row["mode"].endswith("min")
                    else slice(steps // 2, None)
                )
                pw, fw = p[window], f[window]
                scale = float(np.mean(np.abs(fw)))
                difference = abs(float(pw.mean() - fw.mean()))
                allowed = max(relative * scale, floor)
                diff = p - f
                result = {
                    "input": number,
                    "mode": row["mode"],
                    "metric": name,
                    "unit": unit,
                    "window": "final_step"
                    if row["mode"].endswith("min")
                    else "last_half",
                    "window_start": steps
                    if row["mode"].endswith("min")
                    else steps // 2 + 1,
                    "window_end": steps,
                    "python_value": float(pw.mean()),
                    "fortran_value": float(fw.mean()),
                    "absolute_difference": difference,
                    "reference_scale": scale,
                    "difference_percent": 100 * difference / scale
                    if scale > floor
                    else "",
                    "reference_below_absolute_floor": scale <= floor,
                    "relative_tolerance": relative,
                    "absolute_floor": floor,
                    "absolute_floor_controls_band": floor >= relative * scale,
                    "allowed_difference": allowed,
                    "within_band": difference <= allowed,
                    "within_1_percent": difference <= max(0.01 * scale, floor)
                    if relative
                    else "",
                    "within_5_percent": difference <= max(0.05 * scale, floor)
                    if relative
                    else "",
                    "within_10_percent": difference <= max(0.10 * scale, floor)
                    if relative
                    else "",
                    "python_full_mean": float(p.mean()),
                    "fortran_full_mean": float(f.mean()),
                    "python_first": p[0],
                    "fortran_first": f[0],
                    "python_final": p[-1],
                    "fortran_final": f[-1],
                    "full_paired_rmse": float(np.sqrt(np.mean(diff**2))),
                    "full_max_absolute_difference": float(np.max(np.abs(diff))),
                    "full_rmse_percent_reference_scale": 100
                    * float(np.sqrt(np.mean(diff**2)))
                    / float(np.mean(np.abs(f)))
                    if np.mean(np.abs(f)) > floor
                    else "",
                    "python_tail_sd": float(p[steps // 2 :].std(ddof=1)),
                    "fortran_tail_sd": float(f[steps // 2 :].std(ddof=1)),
                    "python_third_quarter_mean": float(
                        p[steps // 2 : 3 * steps // 4].mean()
                    ),
                    "python_fourth_quarter_mean": float(p[3 * steps // 4 :].mean()),
                    "fortran_third_quarter_mean": float(
                        f[steps // 2 : 3 * steps // 4].mean()
                    ),
                    "fortran_fourth_quarter_mean": float(f[3 * steps // 4 :].mean()),
                }
                metrics.append(result)
            if row["mode"] == "metric-dyn":
                for implementation, energy in zip(
                    ("python", "fortran"), plot_data["total_energy"]
                ):
                    conservation.append(
                        {
                            "input": number,
                            "implementation": implementation,
                            "first_total_energy_Ry": energy[0],
                            "final_total_energy_Ry": energy[-1],
                            "relative_energy_range_percent": 100
                            * float(np.ptp(energy))
                            / float(np.mean(np.abs(energy))),
                            "relative_endpoint_drift_percent": 100
                            * float(energy[-1] - energy[0])
                            / abs(float(energy[0])),
                        }
                    )
            write_csv(reference / "observables.csv", normalized)
            completed.append(number)
            fig, axes = plt.subplots(3, 2, figsize=(11, 10), constrained_layout=True)
            for ax, name in zip(
                axes.flat,
                [
                    "atomic_potential_energy",
                    "total_energy",
                    "volume",
                    "pressure",
                    "running_extended_temperature",
                    "cell_kinetic_energy",
                ],
            ):
                p, f = plot_data[name]
                ax.plot(time, p, label="Python", linewidth=0.85)
                ax.plot(time, f, label="Original Fortran", linewidth=0.85, alpha=0.8)
                ax.set(
                    title=name.replace("_", " "),
                    xlabel="Time [rydberg_time]",
                    ylabel=METRICS[name][3],
                )
                ax.axvspan(time[steps // 2], time[-1], color="grey", alpha=0.08)
                ax.ticklabel_format(style="sci", axis="y", scilimits=(-3, 4))
                ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
                ax.grid(alpha=0.2)
            axes[0, 0].legend()
            fig.suptitle(
                f"Input {number}: {row['mode']}, {particles} particles, all {steps} steps"
            )
            pdf.savefig(fig)
            plt.close(fig)
    write_csv(output / "metrics.csv", metrics)
    write_csv(output / "consistency_checks.csv", checks)
    write_csv(output / "energy_conservation.csv", conservation)
    lines = [
        "# Original celq.f comparison results",
        "",
        f"Completed reference workloads: {len(completed)}/8.",
        "",
        "All comparisons below use the predeclared engineering bands and absolute floors. They do not establish statistical equivalence or convergence. These include correlated/invariant observables, not independent tests of equivalence. Dynamic comparisons use last-half means; minimizations use the final step. Percent differences use mean absolute reference magnitude in that window, and are suppressed when the reference is below its absolute floor.",
        "",
        "The maximum absolute difference and paired RMSE over each full trajectory are in metrics.csv; diagnostics.pdf plots all steps. Each pair uses step * input dt after the move/correct operation and before controller events. tv is running extended temperature, not instantaneous atomic temperature.",
        "",
        "Fortran car has only a header; it does not supply a particle trajectory for position-by-position comparisons. io is retained as a historical restart, with its documented incomplete/damaged fields; it is not used as a native checkpoint.",
        "",
        "| Input | Mode | Metrics inside declared bands | Outside bands |",
        "| --- | --- | --- | --- |",
    ]
    for row in read_csv(output / "run_summary.csv"):
        subset = [m for m in metrics if m["input"] == int(row["input"])]
        failed = [m["metric"] for m in subset if not m["within_band"]]
        lines.append(
            f"| {row['input']} | {row['mode']} | {sum(m['within_band'] for m in subset)}/{len(subset)} | {', '.join(failed) or ('none' if subset else row['status'])} |"
        )
    floor_passes = [
        m
        for m in metrics
        if m["within_band"]
        and m["absolute_difference"] > m["relative_tolerance"] * m["reference_scale"]
        and m["relative_tolerance"] > 0
    ]
    lines.extend(
        [
            "",
            "## Absolute-floor comparisons",
            "",
            "These rows pass through the declared absolute floor while exceeding the relative band. Near-zero references make percentages unstable.",
            "",
            "| Input | Observable | Absolute difference | Allowed floor | Unit |",
            "| --- | --- | ---: | ---: | --- |",
        ]
    )
    for m in floor_passes:
        lines.append(
            f"| {m['input']} | {m['metric']} | {m['absolute_difference']:.8g} | {m['absolute_floor']:.8g} | {m['unit']} |"
        )
    structure_path = output / "final_structures.csv"
    if structure_path.exists():
        lines.extend(
            [
                "",
                "## Independent final-structure comparison",
                "",
                "Computed from full-run io geometry and native checkpoints; relative fractional positions are origin-aligned and wrapped. No trajectory was integrated. Single-atom internal structure is uninformative. Thermal trajectories are excluded from atom-by-atom comparison.",
                "",
                "| Input | Cell matrix difference % | Max length difference [bohr] | Relative fractional position error |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        for row in read_csv(structure_path):
            pos = row["position_max_minimum_image_error"]
            position = f"{float(pos):.6g}" if pos else "single atom"
            lines.append(
                f"| {row['input']} | {100 * float(row['cell_frobenius_relative_error']):.6g} | {float(row['cell_length_max_error [bohr]']):.6g} | {position} |"
            )
    lines.extend(
        [
            "",
            "## Energy drift in the original dynamics workloads",
            "",
            "The original timesteps are preserved. Drift agreement does not establish energy conservation. Quenched and thermostatted runs are omitted from this conservation comparison.",
            "",
            "| Input | Implementation | Energy range / mean absolute energy % | Endpoint drift / absolute initial energy % |",
            "| --- | --- | ---: | ---: |",
        ]
    )
    for row in conservation:
        lines.append(
            f"| {row['input']} | {row['implementation']} | {row['relative_energy_range_percent']:.6g} | {row['relative_endpoint_drift_percent']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Detailed primary comparisons",
            "",
            "| Input | Observable | Window | Python | Fortran | Difference % | Inside band |",
            "| --- | --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for m in metrics:
        percent = (
            f"{m['difference_percent']:.6g}"
            if m["difference_percent"] != ""
            else "below absolute floor"
        )
        lines.append(
            f"| {m['input']} | {m['metric']} [{m['unit']}] | {m['window']} | {m['python_value']:.9g} | {m['fortran_value']:.9g} | {percent} | {m['within_band']} |"
        )
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")
    with (directory / "README.md").open("a") as stream:
        stream.write(
            f"\nFinal status: {len(completed)}/8 complete reference workloads; {sum(m['within_band'] for m in metrics)}/{len(metrics)} selected summary checks within their bands or absolute floors; {sum(c['passed'] for c in checks)}/{len(checks)} printed-tape consistency checks passed. See outputs/REPORT.md for limitations and the full numerical results.\n"
        )
    hashes(output, output / "output_hashes.csv")
    print(
        f"Analysis: {len(completed)}/8 complete; {sum(m['within_band'] for m in metrics)}/{len(metrics)} metric comparisons inside declared bands",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-directory", type=Path, default=Path("../vcsmd"))
    parser.add_argument("--python-suite", type=Path)
    parser.add_argument("--examples", type=Path, default=Path("examples/legacy"))
    parser.add_argument("--output-root", type=Path, default=Path("runs"))
    parser.add_argument("--analyze", type=Path)
    args = parser.parse_args()
    if args.analyze:
        analyze(args.analyze)
        return 0
    if args.python_suite is None:
        parser.error("--python-suite is required when executing the reference")
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise RuntimeError("gfortran is required")
    stamp = datetime.now().astimezone().strftime("%Y-%m-%d-%H")
    args.output_root.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while True:
        directory = (
            args.output_root
            / f"vcsmd-comparison-celq-eight-inputs-results-{stamp}-attempt-{attempt:02d}"
        )
        try:
            directory.mkdir()
            break
        except FileExistsError:
            attempt += 1
    inputs, output = directory / "inputs", directory / "outputs"
    inputs.mkdir()
    output.mkdir()
    (directory / "README.md").write_text(
        PROTOCOL
        + f"\nCreated: {stamp}. Attempt: {attempt}.\n\nPython provenance: copied from the existing suite {args.python_suite.name}.\n"
    )
    (inputs / "fortran_source").mkdir()
    for name in ("celq.f", "sip"):
        shutil.copyfile(args.legacy_directory / name, inputs / "fortran_source" / name)
    shutil.copytree(args.examples, inputs / "legacy_inputs")
    shutil.copyfile(__file__, inputs / "compare_fortran.py")
    (inputs / "python_runs").mkdir()
    original = read_csv(args.python_suite / "outputs/summary.csv")
    selected = [r for r in original if r["format"] == "yaml"]
    if sorted(int(r["input"]) for r in selected) != list(range(1, 9)):
        raise ValueError("Eight saved Python YAML runs required")
    for row in selected:
        n = int(row["input"])
        source = args.python_suite / row["run"]
        config = load_config(source / "inputs/configuration.json")
        legacy = config_from_mapping(
            parse_legacy_input(
                (inputs / "legacy_inputs" / f"input-{n:02d}.inp").read_text()
            )
        )
        if config_to_mapping(config) != config_to_mapping(legacy):
            raise ValueError(f"Input {n}: Python and legacy configurations differ")
        if row["status"] != "complete" or int(row["completed_steps"]) != config.steps:
            raise ValueError(f"Input {n}: prior Python run incomplete")
        shutil.copytree(source, inputs / "python_runs" / f"input-{n:02d}")
    hashes(inputs, inputs / "source_hashes.csv")
    version = subprocess.run(
        [compiler, "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    (output / "build").mkdir()
    # Keep original fixed-form line length; do not activate code beyond column 72.
    command = [
        compiler,
        "-std=legacy",
        "-O2",
        "-o",
        "celq",
        "../../inputs/fortran_source/celq.f",
    ]
    built = subprocess.run(
        command, cwd=output / "build", capture_output=True, text=True, check=False
    )
    (output / "build/compiler.log").write_text(
        version
        + "\nCommand: gfortran "
        + " ".join(command[1:])
        + "\n"
        + built.stdout
        + built.stderr
    )
    if built.returncode:
        print("Fortran compilation failed; see outputs/build/compiler.log", flush=True)
        return 1
    print(f"Compiled unmodified celq.f: {version}; suite {directory.name}", flush=True)
    (output / "fortran").mkdir()
    summaries = []
    for row in selected:
        n = int(row["input"])
        path = (
            output / "fortran" / f"vcsmd-celq-input-{n:02d}-results-{stamp}-attempt-01"
        )
        path.mkdir()
        for name, source in [
            ("inp", inputs / "legacy_inputs" / f"input-{n:02d}.inp"),
            ("sip", inputs / "fortran_source/sip"),
        ]:
            (path / name).symlink_to(os.path.relpath(source, path))
        (path / "README.md").write_text(
            f"# Full Fortran input {n}\n\nCreated {stamp}, attempt 1. Executes the original {row['particles']}-particle, {row['requested_steps']}-step workload. Inputs inp and sip are symlinks to copied suite input snapshots. Output tapes e/eal/ave/p/avec/tv, out, car, io and execution.log are preserved. Binary and compiler flags are in ../../build. No source or simulation parameter is modified.\n"
        )
        print(
            f"Fortran input {n}: starting all {row['requested_steps']} steps",
            flush=True,
        )
        with (path / "execution.log").open("w") as stream:
            execution = subprocess.run(
                [str((output / "build/celq").resolve())],
                cwd=path,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        status, error = "complete", ""
        try:
            if execution.returncode:
                raise ValueError(f"process exit {execution.returncode}")
            for name, count in {"e": 4, "eal": 6, "ave": 2, "p": 3, "tv": 2}.items():
                tape(path / name, count, int(row["requested_steps"]))
            if row["mode"] != "fixed-dyn":
                tape(path / "avec", 6, int(row["requested_steps"]))
        except (ValueError, OSError, IndexError) as exc:
            status, error = "failed", str(exc)
        summaries.append(
            {
                "input": n,
                "mode": row["mode"],
                "particles": row["particles"],
                "steps": row["requested_steps"],
                "status": status,
                "error": error,
                "fortran_directory": path.relative_to(output).as_posix(),
            }
        )
        write_csv(output / "run_summary.csv", summaries)
        print(
            f"Fortran input {n}: {status}" + (f" ({error})" if error else ""),
            flush=True,
        )
    analyze(directory)
    return 0 if all(r["status"] == "complete" for r in summaries) else 1


if __name__ == "__main__":
    raise SystemExit(main())
