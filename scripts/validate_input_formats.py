"""Run all eight shipped YAML and CLI-imported legacy workloads in full.

No reduced trajectories are created. Run from the repository root with
``uv run python scripts/validate_input_formats.py``.
"""

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

from vcsmd import RunStatus, config_from_mapping, config_to_mapping, load_config, run
from vcsmd.compat import parse_legacy_input, split_legacy_examples


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def checkpoints_equal(first: Path, second: Path) -> bool:
    with (
        np.load(first, allow_pickle=False) as a,
        np.load(second, allow_pickle=False) as b,
    ):
        return set(a.files) == set(b.files) and all(
            np.array_equal(a[key], b[key]) for key in a.files
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", type=Path, default=Path("examples"))
    parser.add_argument("--output-root", type=Path, default=Path("runs"))
    args = parser.parse_args()
    stamp = datetime.now().astimezone().strftime("%Y-%m-%d-%H")
    args.output_root.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while True:
        directory = args.output_root / (
            f"vcsmd-validation-yaml-legacy-inputs-results-{stamp}-attempt-{attempt:02d}"
        )
        try:
            directory.mkdir()
            break
        except FileExistsError:
            attempt += 1
    inputs, outputs = directory / "inputs", directory / "outputs"
    inputs.mkdir()
    outputs.mkdir()
    copied = inputs / "examples"
    shutil.copytree(args.examples, copied)
    shutil.copyfile(__file__, inputs / "validate_input_formats.py")
    write_csv(
        inputs / "source_hashes.csv",
        [
            {
                "file": p.relative_to(inputs).as_posix(),
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
            for p in sorted(inputs.rglob("*"))
            if p.is_file()
        ],
    )
    (directory / "README.md").write_text(
        "# Full YAML and legacy input validation\n\n"
        f"Created: {stamp}. Purpose: execute all eight shipped YAML inputs and all "
        "eight CLI-imported Fortran-format inputs at their original particle counts, "
        "parameters and total step counts (16 full runs, 18,400 integration steps).\n\n"
        "Inputs: examples/ is a copied snapshot of the supplied example directory; "
        "source_hashes.csv records its provenance together with this validation script. "
        "imported/ contains native YAML produced by the actual import-legacy CLI. "
        "Each child run copies its input and executed Python source and records "
        "dependency versions and hashes. Physical input parameters are preserved. "
        "Historical input records contain no random seed; the importer supplies "
        "native seed 119, matching the shipped YAML files.\n\n"
        "Outputs: summary.csv records completion; checks.csv records equivalence, "
        "finite observables and exact paired checkpoint/CSV comparisons. import-*.log "
        "records CLI conversion output. runs/ contains complete native results, CSV "
        "streams, checkpoints, provenance and status. All outputs are retained. "
        "No temporary or reduced simulation cases are created.\n\n"
        "Reproduce from the repository root: uv run python scripts/validate_input_formats.py. "
        "Use --examples to select the preserved inputs/examples snapshot and "
        "--output-root to choose a new output parent. This validates the Python "
        "legacy adapter, not numerical parity with a Fortran executable.\n\n"
        "Completion is confirmed only by the final status appended below after all "
        "checks pass. Otherwise the suite is incomplete; consult summary.csv and "
        "checks.csv for recorded progress and failures.\n"
    )
    (inputs / "imported").mkdir()
    blocks = split_legacy_examples((copied / "original-inputs.txt").read_text())
    if set(blocks) != set(range(1, 9)):
        raise ValueError("All eight original input blocks are required")
    checks: list[dict] = []
    rows: list[dict] = []

    def check(name: str, passed: bool) -> None:
        checks.append({"check": name, "passed": bool(passed)})
        write_csv(outputs / "checks.csv", checks)
        if not passed:
            raise ValueError(f"Validation failed: {name}")

    for number in range(1, 9):
        stem = f"input-{number:02d}"
        legacy = copied / "legacy" / f"{stem}.inp"
        check(
            f"Input {number} original legacy records preserved",
            legacy.read_text() == blocks[number].lstrip("\r\n"),
        )
        imported = inputs / "imported" / f"{stem}.yaml"
        converted = subprocess.run(
            [
                sys.executable,
                "-m",
                "vcsmd",
                "import-legacy",
                str(legacy),
                str(imported),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        (outputs / f"import-{number:02d}.log").write_text(
            converted.stdout + converted.stderr
        )
        check(f"Input {number} CLI import succeeds", converted.returncode == 0)
        paths = {"yaml": copied / f"{stem}.yaml", "legacy": imported}
        reference = config_to_mapping(
            config_from_mapping(parse_legacy_input(blocks[number]))
        )
        for label, path in {"toml": copied / f"{stem}.toml", **paths}.items():
            check(
                f"Input {number} {label} matches original configuration",
                config_to_mapping(load_config(path)) == reference,
            )
        if number == 1:
            check(
                "Input 1 JSON matches original configuration",
                config_to_mapping(load_config(copied / f"{stem}.json")) == reference,
            )
        reports = {}
        for label, path in paths.items():
            config = load_config(path)
            print(f"Input {number} {label}: {config.steps} original steps", flush=True)

            def progress(
                completed: int, total: int, number: int = number, label: str = label
            ) -> None:
                if completed % 500 == 0 or completed == total:
                    print(f"Input {number} {label}: {completed}/{total}", flush=True)

            report = run(
                config,
                source_file=path,
                output_root=outputs / "runs",
                purpose=f"{stem}-{label}",
                checkpoint_interval=100,
                progress=progress,
            )
            reports[label] = report
            rows.append(
                {
                    "input": number,
                    "format": label,
                    "mode": config.mode.value,
                    "particles": len(config.structure.species),
                    "requested_steps": config.steps,
                    "completed_steps": report.completed_steps,
                    "status": report.status.value,
                    "error": report.error or "",
                    "run": report.directory.relative_to(directory).as_posix(),
                }
            )
            write_csv(outputs / "summary.csv", rows)
            check(
                f"Input {number} {label} full workload complete",
                report.status is RunStatus.COMPLETE
                and report.completed_steps == config.steps,
            )
            with (report.directory / "outputs" / "observables.csv").open() as stream:
                records = list(csv.DictReader(stream))
            check(
                f"Input {number} {label} finite observables at every step",
                len(records) == config.steps
                and all(
                    np.isfinite(float(value))
                    for row in records
                    for value in row.values()
                ),
            )
        first = reports["yaml"].directory / "outputs"
        second = reports["legacy"].directory / "outputs"
        check(
            f"Input {number} exact final checkpoint equality",
            checkpoints_equal(first / "checkpoint.npz", second / "checkpoint.npz"),
        )
        first_csv = {p.name for p in first.glob("*.csv")}
        second_csv = {p.name for p in second.glob("*.csv")}
        check(
            f"Input {number} exact complete CSV stream equality",
            bool(first_csv)
            and first_csv == second_csv
            and all(
                (first / name).read_bytes() == (second / name).read_bytes()
                for name in first_csv
            ),
        )
    with (directory / "README.md").open("a") as stream:
        stream.write(
            f"\nFinal status: all {len(rows)} full runs and {len(checks)} checks passed.\n"
        )
    print(
        f"Validation complete: {len(rows)} full runs, {len(checks)} checks; {directory.name}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
