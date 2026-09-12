"""One-way readers for the historical VCSMD input and output formats.

The modules in :mod:`vcsmd.compat` deliberately have no dependency on the
numerical implementation.  They turn the fixed-format files written by the
Fortran program into the native, JSON/CSV-friendly data contract used by the
application layer.  Historical calculation codes and filenames are confined
to this package.
"""

from .legacy import (
    ConversionReport,
    LegacyFormatError,
    convert_legacy_run,
    import_legacy_input,
    parse_legacy_input,
    split_legacy_examples,
)

__all__ = [
    "ConversionReport",
    "LegacyFormatError",
    "convert_legacy_run",
    "import_legacy_input",
    "parse_legacy_input",
    "split_legacy_examples",
]
