"""Check installed APIs with supplied configurations and completed full-run state.

This performs no time integration and creates no simulation workload. Temporary
serialization files are removed when their checks finish.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import importlib.metadata
import pkgutil
import platform
import subprocess
import sys
import tempfile
from dataclasses import fields, is_dataclass
from pathlib import Path

import numpy as np

import vcsmd
from vcsmd.cli import main as cli_main
from vcsmd.compat import parse_legacy_input, split_legacy_examples
from vcsmd.config import prepare
from vcsmd.io import (
    config_from_mapping,
    load_checkpoint,
    load_config,
    save_checkpoint,
)


def equal(first: object, second: object) -> bool:
    if type(first) is not type(second):
        return False
    if isinstance(first, np.ndarray):
        return np.array_equal(first, second)
    if is_dataclass(first):
        return all(
            equal(getattr(first, item.name), getattr(second, item.name))
            for item in fields(first)
        )
    if isinstance(first, tuple):
        return len(first) == len(second) and all(
            equal(a, b) for a, b in zip(first, second)
        )
    return first == second


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("validation_run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []

    def record(check: str, passed: bool, detail: str = "") -> None:
        rows.append({"check": check, "passed": bool(passed), "detail": detail})
        print(f"{check}: {'PASS' if passed else 'FAIL'} {detail}")

    versions = {"python": platform.python_version()}
    for name in ("vcsmd", "numpy", "pint", "PyYAML", "tomli-w"):
        versions[name] = importlib.metadata.version(name)
    record("Installed dependency versions", True, repr(versions))
    modules = [
        module.name
        for module in pkgutil.walk_packages(vcsmd.__path__, vcsmd.__name__ + ".")
        if not module.name.endswith(".__main__")
    ]
    for name in modules:
        importlib.import_module(name)
    record("All installed package modules import", True, f"{len(modules)} modules")
    result = subprocess.run(
        [sys.executable, "-m", "vcsmd", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    record("Installed module CLI", result.returncode == 0 and "resume" in result.stdout)

    prepared = [
        prepare(
            load_config(args.validation_run / "inputs" / f"format-equivalence.{ext}")
        )
        for ext in ("json", "yaml", "toml")
    ]
    record(
        "JSON/YAML/TOML normalized model equivalence",
        all(equal(prepared[0], x) for x in prepared[1:]),
    )
    with (args.validation_run / "outputs" / "summary.csv").open(newline="") as stream:
        second = next(row for row in csv.DictReader(stream) if row["input"] == "2")
    checkpoint = args.validation_run / second["run"] / "outputs" / "checkpoint.npz"
    model, state = load_checkpoint(checkpoint)
    record(
        "Public checkpoint returns dimensional state",
        state.cell.check("[length]") and state.time.check("[time]"),
    )
    record(
        "Public state arrays are read-only",
        not state.cell.magnitude.flags.writeable
        and not state.fractional_positions.flags.writeable,
    )

    with tempfile.TemporaryDirectory(prefix="vcsmd-serialization-") as temporary:
        directory = Path(temporary)
        copy = directory / "checkpoint.npz"
        save_checkpoint(copy, model, state)
        restored_model, restored_state = load_checkpoint(copy)
        record(
            "Public checkpoint roundtrip preserves all fields",
            equal(model, restored_model)
            and equal(state.numerical, restored_state.numerical),
        )
        record(
            "Loaded state owns independent arrays",
            not np.shares_memory(state.cell.magnitude, restored_state.cell.magnitude),
        )
        original = split_legacy_examples(
            (args.validation_run / "inputs" / "original-inputs.txt").read_text()
        )[1]
        source = directory / "original-input.txt"
        destination = directory / "configuration.toml"
        source.write_text(original)
        exit_code = cli_main(["import-legacy", str(source), str(destination)])
        record(
            "CLI imports complete original input",
            exit_code == 0
            and equal(
                prepare(load_config(destination)),
                prepare(config_from_mapping(parse_legacy_input(original))),
            ),
        )

    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["check", "passed", "detail"])
        writer.writeheader()
        writer.writerows(rows)
    return 0 if all(row["passed"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
