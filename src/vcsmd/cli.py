"""Command-line boundary for native execution and isolated legacy conversion."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from .execution import RunStatus, resume, run
from .io import config_from_mapping, load_config, save_config


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vcsmd", description="Unit-aware variable-cell molecular dynamics"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    native = subparsers.add_parser(
        "run", help="Run a JSON, YAML, or TOML configuration"
    )
    native.add_argument("configuration", type=Path)
    native.add_argument("--output-root", type=Path, default=Path("runs"))
    native.add_argument("--checkpoint-interval", type=int, default=100)
    continuation = subparsers.add_parser(
        "resume", help="Continue an exact native checkpoint"
    )
    continuation.add_argument("checkpoint", type=Path)
    continuation.add_argument(
        "--steps", type=int, required=True, help="Additional integration steps"
    )
    continuation.add_argument("--output-root", type=Path, default=Path("runs"))
    continuation.add_argument("--output-interval", type=int, default=1)
    continuation.add_argument("--checkpoint-interval", type=int, default=100)
    converter = subparsers.add_parser(
        "convert-legacy", help="Convert an old run directory into native datasets"
    )
    converter.add_argument("source", type=Path)
    converter.add_argument("--destination", type=Path)
    importer = subparsers.add_parser(
        "import-legacy", help="Convert an old input into a native configuration"
    )
    importer.add_argument("source", type=Path)
    importer.add_argument(
        "destination", type=Path, help="New .json, .yaml, or .toml file"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            report = run(
                load_config(args.configuration),
                output_root=args.output_root,
                source_file=args.configuration,
                checkpoint_interval=args.checkpoint_interval,
            )
        elif args.command == "resume":
            report = resume(
                args.checkpoint,
                steps=args.steps,
                output_root=args.output_root,
                output_interval=args.output_interval,
                checkpoint_interval=args.checkpoint_interval,
            )
        elif args.command == "import-legacy":
            from .compat import import_legacy_input

            if args.destination.exists():
                raise FileExistsError("Destination configuration already exists")
            save_config(
                config_from_mapping(import_legacy_input(args.source)), args.destination
            )
            print("Native configuration written.")
            return 0
        else:
            from .compat import convert_legacy_run

            destination = args.destination
            if destination is None:
                stamp = datetime.now().astimezone().strftime("%Y-%m-%d-%H")
                for attempt in range(1, 10000):
                    candidate = (
                        Path("runs")
                        / f"vcsmd-conversion-legacy-native-results-{stamp}-attempt-{attempt:02d}"
                    )
                    if not candidate.exists():
                        destination = candidate
                        break
            converted = convert_legacy_run(args.source, destination)
            print(
                f"Conversion completed: {len(converted.problems)} reported limitations or malformed records."
            )
            return 0
        print(
            f"{report.status.value}: {report.completed_steps}/{report.requested_steps} steps; {report.directory.name}"
        )
        if report.error:
            print(report.error)
        return 0 if report.status is RunStatus.COMPLETE else 1
    except (ValueError, TypeError, OSError) as exc:
        # File exceptions can embed identifying absolute paths; expose their
        # category only. Validation errors use field names and physical values.
        detail = (
            type(exc).__name__
            if isinstance(exc, OSError)
            else str(exc).replace(str(Path.home()), "~")
        )
        parser.exit(2, f"vcsmd: {detail}\n")
    return 2
