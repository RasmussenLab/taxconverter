# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Taxconverter is a Python CLI tool that converts outputs from different taxonomic classifiers (MMseqs2, Centrifuge, Kraken2, Metabuli, MetaMaps) into a unified format. The tool standardizes taxonomic annotations to a common `contigs\tpredictions` format, where predictions are full lineages from domain to species concatenated with semicolons.

The package uses NCBI taxonomy IDs (for Centrifuge/Kraken2/MetaMaps) or GTDB identifiers (for Metabuli/MMseqs2) and resolves complete lineages using a reference file (`data/clades.tsv`).

IMPORTANT: THIS PROJECT IS UNDERGOING A REFACTOR. Hence, the observed patterns may not match the ideal patterns written in this Claude.md.

## Key Architecture

### One module per source format
This package should have one module per source, e.g. an `mmseqs` module, a `kraken` module etc.
However, this may not be implemented yet, and code from several sources may be in the same module.

The CLI lives in `taxconverter/__main__.py`

### CLI Structure

The tool uses argparse with subcommands for each classifier:
- `taxconverter centrifuge` - Converts Centrifuge output
- `taxconverter kraken2` - Converts Kraken2 output
- `taxconverter metabuli` - Converts Metabuli output (requires 2 files: classifications and report)
- `taxconverter metamaps` - Converts MetaMaps output
- `taxconverter mmseqs2` - Converts MMseqs2 output

### Core Data Flow
The ideal core data flows is listed below. Note that this may not be what is actually implemented, but rather what the refactor works towards:

1. Parse the input data to a common format, the `Clades` object.
1a. Optionally parse the NCBI lineage in order to resolve the linege from NCBI names/identifiers
2. Write the output in a unified format, either in MMseqs format, or in an internal "Vamb" format (which is a TSV file)

### Metabuli Special Handling

Metabuli requires custom parsing because:
- It uses dataset-specific GTDB identifiers rather than NCBI tax IDs
- The report file format changed across versions (header was added in commit a17debb)
- The lineage must be reconstructed from indented clade names in the report file
- The parser (`metabuli_from_iter()`) is strict about format validation to catch version incompatibilities

The algorithm tracks canonical ranks by indentation depth and builds lineages incrementally as it reads the report.

## Development Commands

### Installation

```bash
# Install in development mode
pip install -e .
```

The package requires Python <3.12 due to pandas 3.0 compatibility constraints.

### Running the Tool

```bash
# General pattern
taxconverter <classifier> -i <input> -o <output>

# Metabuli (requires two files)
taxconverter metabuli -c classifications.tsv -r report.tsv -o output.tsv

# With MMseqs2 format (for older Taxometer compatibility)
taxconverter kraken2 -i input.tsv -o output.tsv --mmseqs-format
```

### Testing

Tests are configured in `pyproject.toml` under `[tool.pytest.ini_options]`. Run with:

```bash
pytest
```

Note: The repository currently has no test files (no `tests/` directory exists).

### Code Quality

The project uses ruff for linting. Configuration in `pyproject.toml` ignores:
- E722 (bare except)
- E501 (line too long)

Run linting with:

```bash
ruff check .
```

## Important Implementation Details

### Version Management

Version is defined as a tuple in `taxconverter/__init__.py`:
```python
__version__ = (1, 1, 0)
```

This is referenced dynamically in `pyproject.toml` via `tool.setuptools.dynamic`.

### Data Files

The NCBI taxonomy file `data/clades.tsv` is packaged with the tool via `[tool.setuptools.package-data]` in `pyproject.toml`. The path is resolved relative to the parent directory in `__main__.py`:

```python
NCBI_LINEAGE = os.path.join(parentdir, 'data/clades.tsv')
```

### Type System
This project ought to use type hints thoroughly, and to avoid untyped containers like pandas dataframes. However, this ideal might not be fully realized yet.

### Error Handling

- Missing tax IDs in NCBI lineage are logged but return empty lineage strings
- Metabuli parsing is strict and raises `ValueError` for format violations
- The decorator pattern ensures all converters apply consistent output formatting
