"""One-way readers for the historical VCSMD input and output formats.

These adapters turn historical fixed-format files into native configuration
mappings and CSV import datasets. Parsing uses no solver operations; conversion
validates generated configurations through the native configuration adapter.
The computational core never imports this package. Historical calculation codes
and filenames are confined to this compatibility boundary.
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
