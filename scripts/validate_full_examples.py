"""Run the eight original workloads without reducing their sizes or durations.

This is an auditable acceptance workflow, not a generator of synthetic cases.
Invoke after installing the project's dev dependencies (for PDF plots).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from dataclasses import fields
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from vcsmd.compat import convert_legacy_run, parse_legacy_input, split_legacy_examples
from vcsmd.execution import RunStatus, resume, run
from vcsmd.io import config_from_mapping, config_to_mapping, load_config, save_config
from vcsmd.io.checkpoint import load_checkpoint


def _json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _table(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    if not rows:
        return {}
    return {key: np.asarray([float(row[key]) for row in rows]) for key in rows[0]}


def _check_state_equal(first: object, second: object) -> bool:
    for item in fields(first):
        a = getattr(first, item.name)
        b = getattr(second, item.name)
        if isinstance(a, np.ndarray):
            if not np.array_equal(a, b):
                return False
        elif a != b:
            return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=Path("examples/original-inputs.txt")
    )
    parser.add_argument("--output-root", type=Path, default=Path("runs"))
    parser.add_argument("--legacy-directory", type=Path)
    args = parser.parse_args()
    stamp = datetime.now().astimezone().strftime("%Y-%m-%d-%H")
    args.output_root.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while True:
        directory = (
            args.output_root
            / f"vcsmd-validation-eight-inputs-results-{stamp}-attempt-{attempt:02d}"
        )
        try:
            directory.mkdir()
            break
        except FileExistsError:
            attempt += 1
    inputs = directory / "inputs"
    output = directory / "outputs"
    inputs.mkdir()
    output.mkdir()
    shutil.copyfile(args.source, inputs / "original-inputs.txt")
    shutil.copyfile(Path(__file__), inputs / "validate_full_examples.py")
    source_text = args.source.read_text()
    examples = split_legacy_examples(source_text)
    if len(examples) != 8:
        raise ValueError("This workflow requires all eight original examples")
    _json(
        inputs / "provenance.json",
        {
            "source": "copied complete original input collection",
            "sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
            "created_at": datetime.now().astimezone().isoformat(timespec="hours"),
            "parameters_modified": False,
        },
    )
    (directory / "README.md").write_text(
        "# Full original-example validation\n\n"
        f"Created: {stamp}. Purpose: execute all eight original examples at their original "
        "particle counts, time steps, temperatures, and total step counts.\n\n"
        "Inputs contain a copied original input collection, its hash/provenance, normalized configurations "
        "and the acceptance script. Each child run also snapshots implementation source and dependencies.\n\n"
        "Outputs contain complete native runs, checkpoint-continuation verification, a CSV summary, "
        "format/purity checks and PDF diagnostics. The continuation workload completes the same original "
        "trajectory from a saved intermediate checkpoint; it is not a shorter replacement workload.\n\n"
        "Existing examples cover fixed-dyn, metric-dyn and metric-min. Other modes receive source/equation "
        "review; no additional trajectories are synthesized. No numerical Fortran baseline is claimed.\n"
    )
    rows = []
    reports = {}
    checks = []
    with PdfPages(output / "diagnostics.pdf") as pdf:
        for number, original in examples.items():
            config = config_from_mapping(parse_legacy_input(original))
            config_path = inputs / f"input-{number:02d}.json"
            save_config(config, config_path)
            before = json.dumps(config_to_mapping(config), sort_keys=True)
            if number == 1:
                equivalent = []
                for suffix in ("json", "yaml", "toml"):
                    filename = inputs / f"format-equivalence.{suffix}"
                    save_config(config, filename)
                    equivalent.append(
                        json.dumps(
                            config_to_mapping(load_config(filename)), sort_keys=True
                        )
                    )
                checks.append(
                    {
                        "check": "JSON/YAML/TOML model equivalence",
                        "passed": len(set(equivalent)) == 1,
                    }
                )
            print(
                f"Input {number}: {config.mode.value}, {len(config.structure.species)} particles, {config.steps} steps",
                flush=True,
            )

            def progress(completed: int, total: int, number: int = number) -> None:
                if completed % 500 == 0 or completed == total:
                    print(f"Input {number}: {completed}/{total}", flush=True)

            report = run(
                config,
                output_root=output / "runs",
                purpose=f"validation-input-{number}",
                source_file=config_path,
                checkpoint_interval=100,
                progress=progress,
            )
            reports[number] = report
            checks.append(
                {
                    "check": f"Input {number} configuration unchanged",
                    "passed": before
                    == json.dumps(config_to_mapping(config), sort_keys=True),
                }
            )
            row = {
                "input": number,
                "mode": config.mode.value,
                "particles": len(config.structure.species),
                "requested_steps": config.steps,
                "completed_steps": report.completed_steps,
                "status": report.status.value,
                "energy_relative_range": "",
                "final_volume_bohr3": "",
                "maximum_atomic_temperature_K": "",
                "error": report.error or "",
                "run": report.directory.relative_to(directory).as_posix(),
            }
            observations = report.directory / "outputs" / "observables.csv"
            data = _table(observations) if observations.exists() else {}
            if data:
                time = data["time [rydberg_time]"]
                energy = data["total_energy [Ry]"]
                row["energy_relative_range"] = float(
                    np.ptp(energy) / max(abs(float(energy[0])), np.finfo(float).tiny)
                )
                row["final_volume_bohr3"] = float(data["volume [bohr^3]"][-1])
                row["maximum_atomic_temperature_K"] = float(
                    np.max(data["atomic_temperature [K]"])
                )
                checks.append(
                    {
                        "check": f"Input {number} finite observables",
                        "passed": all(np.all(np.isfinite(v)) for v in data.values()),
                    }
                )
                fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
                axes[0, 0].plot(time, energy)
                axes[0, 0].set_ylabel("Extended energy [Ry]")
                axes[0, 1].plot(time, data["volume [bohr^3]"])
                axes[0, 1].set_ylabel("Volume [bohr³]")
                axes[1, 0].plot(time, data["atomic_temperature [K]"], label="Atomic")
                axes[1, 0].plot(
                    time,
                    data["running_extended_temperature [K]"],
                    label="Running extended",
                )
                axes[1, 0].set_ylabel("Temperature [K]")
                axes[1, 0].legend()
                axes[1, 1].plot(time, data["pressure [Ry/bohr^3]"])
                axes[1, 1].set_ylabel("Pressure [Ry/bohr³]")
                for ax in axes.flat:
                    ax.set_xlabel("Time [Rydberg time]")
                fig.suptitle(
                    f"Original input {number}: {config.mode.value}; {report.status.value}"
                )
                pdf.savefig(fig)
                plt.close(fig)
            rows.append(row)
            with (output / "summary.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
            print(
                f"Input {number}: {report.status.value}"
                + (f" ({report.error})" if report.error else ""),
                flush=True,
            )
    # Complete the ORIGINAL input 2 workload from its halfway checkpoint and
    # compare its final state with the uninterrupted, full-length execution.
    reference = reports[2]
    if reference.status is RunStatus.COMPLETE:
        middle = reference.directory / "outputs" / "checkpoints" / "step-00001000.npz"
        continuation = resume(
            middle,
            steps=1000,
            output_root=output / "runs",
            purpose="validation-input-2-full-continuation",
            checkpoint_interval=100,
        )
        exact = False
        if continuation.status is RunStatus.COMPLETE:
            _, expected = load_checkpoint(
                reference.directory / "outputs" / "checkpoint.npz"
            )
            _, actual = load_checkpoint(
                continuation.directory / "outputs" / "checkpoint.npz"
            )
            exact = _check_state_equal(expected, actual)
        checks.append(
            {
                "check": "Input 2 full-length checkpoint continuation exactly equal",
                "passed": exact,
            }
        )
    else:
        checks.append(
            {
                "check": "Input 2 full-length checkpoint continuation unavailable",
                "passed": False,
            }
        )
    if args.legacy_directory is not None:
        converted = convert_legacy_run(
            args.legacy_directory,
            output / f"vcsmd-conversion-original-repository-results-{stamp}-attempt-01",
        )
        _json(output / "legacy-conversion-availability.json", converted.as_dict())
    with (output / "checks.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(checks)
    completed = sum(row["status"] == "complete" for row in rows)
    _json(
        output / "validation_metadata.json",
        {
            "full_runs_completed": completed,
            "full_runs_requested": 8,
            "checks_passed": all(check["passed"] for check in checks),
            "executed_modes": sorted({row["mode"] for row in rows}),
            "fortran_numerical_baseline": "unavailable",
            "synthetic_or_reduced_simulations_created": False,
        },
    )
    print(
        f"Validation complete: {completed}/8 complete; artifacts: {directory.name}",
        flush=True,
    )
    return 0 if completed == 8 and all(check["passed"] for check in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
