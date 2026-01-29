# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

Taxconverter is a Python CLI tool that converts outputs from 5 different taxonomic classifiers into a unified format:
- **Centrifuge** (NCBI taxonomy)
- **Kraken2** (NCBI taxonomy)
- **MetaMaps** (NCBI taxonomy)
- **Metabuli** (GTDB identifiers)
- **MMseqs2** (GTDB identifiers)

The tool standardizes taxonomic annotations to a common TSV format with full lineages from domain to species.

## Architecture

### Module Structure (One Module Per Parser)

The package has one module per source format:
- `taxconverter/centrifuge.py` - Centrifuge parser
- `taxconverter/kraken.py` - Kraken2 parser
- `taxconverter/metabuli.py` - Metabuli parser
- `taxconverter/metamaps.py` - MetaMaps parser
- `taxconverter/mmseqs.py` - MMseqs2 parser
- `taxconverter/common.py` - Core types and NCBIRanks
- `taxconverter/__main__.py` - CLI entry point

### Core Type System

The package uses a strictly-typed dataclass hierarchy defined in `taxconverter/common.py`:

```python
# Intermediate representation (parser output)
@dataclass
class UnvalidatedIntAnnotation:
    contig_name: str
    leaf_clade: Optional[int]  # None if unclassified

# NCBI-based annotation (after resolving lineage)
@dataclass
class NCBIAnnotation:
    contig_name: str
    clades: list[NCBIID]  # Ordered list from domain to species

# Final generic annotation (output format)
@dataclass
class GenericAnnotation:
    contig_name: str
    clades: list[NCBIIdentifier]  # Ordered lineage as strings

    def to_string(self) -> str:
        # Returns "contig_name\tlineage" format

    def get_mmseqs_rank(self) -> str:
        # Returns rank based on lineage length

# Helper types
NCBIID = NewType("NCBIID", int)
NCBIIdentifier = str wrapper that validates (no semicolons)
Rank = NewType("Rank", int)  # 0-6 for canonical ranks
```

**Canonical ranks** (7 total):
```python
["domain", "phylum", "class", "order", "family", "genus", "species"]
```

### Data Flow

**NCBI-based parsers** (Centrifuge, Kraken2, MetaMaps):
1. Parse file → `list[UnvalidatedIntAnnotation]` (contig + taxid)
2. Call `ncbi.annotation_from_int()` → `list[NCBIAnnotation]` (resolved NCBI IDs)
3. Call `ncbi.generic_annotation()` → `list[GenericAnnotation]` (names as strings)
4. Write output

**Direct lineage parsers** (Metabuli, MMseqs2):
1. Parse file → `list[GenericAnnotation]` (already have lineage names)
2. Write output

### NCBIRanks Class

Located in `taxconverter/common.py`, this class maps NCBI taxonomy IDs to full lineages:

```python
class NCBIRanks:
    child_data: dict[NCBIID, tuple[NCBIID, Optional[Rank], NCBIIdentifier]]
    # Maps: child_id → (parent_id, rank, name)

    @classmethod
    def from_file(cls, path: Path) -> Self:
        # Loads from data/clades.tsv or data/clades.tsv.gz

    def annotation_from_int(self, int_annotation: UnvalidatedIntAnnotation) -> NCBIAnnotation:
        # Resolves taxid to full lineage by walking up parent chain

    def generic_annotation(self, ncbi_annotation: NCBIAnnotation) -> GenericAnnotation:
        # Converts NCBI IDs to names
```

**Data file**: `data/clades.tsv.gz` (or `.tsv` uncompressed)
- Packaged with the tool via `pyproject.toml` setuptools config
- Format: `child_id\tchild_rank\tparent_id\tname`
- Generated from https://github.com/RasmussenLab/misc_scripts/blob/master/parse_ncbi_tax.jl
- Handles both canonical ranks and intermediate ranks (e.g., infraorder)

## Parser Functions

### Centrifuge (`taxconverter/centrifuge.py`)

```python
def parse_centrifuge(path: Path, lines: Iterator[str], ncbi: NCBIRanks) -> list[NCBIAnnotation]
```

**Input format**: 8 tab-separated columns with header
```
readID	seqID	taxID	score	2ndBestScore	hitLength	queryLength	numMatches
```

**Behavior**:
- Validates header exactly
- taxID=0 means unclassified (→ empty clades list)
- Resolves taxIDs via NCBIRanks

### Kraken2 (`taxconverter/kraken.py`)

```python
def parse_kraken(path: Path, lines: Iterable[str], ncbi: NCBIRanks) -> list[NCBIAnnotation]
```

**Input format**: 5 tab-separated columns, **no header**
```
C/U	readID	taxID	length	lca_mapping
```

**Behavior**:
- C = classified, U = unclassified
- No header validation (starts parsing immediately)
- taxID=0 means unclassified

### Metabuli (`taxconverter/metabuli.py`)

```python
def parse_metabuli_files(classification_path: Path, report_path: Path) -> list[GenericAnnotation]
```

**Requires TWO files**:
1. **Classifications**: 7 columns (optional header as of commit a17debb)
   ```
   db_id	query_id	ref_id	ref_len	score	rank	counts
   ```
   Column 1: `is_classified` (0 or 1)

2. **Report**: 6 columns, indented hierarchy
   ```
   percentage	reads	direct_reads	rank	taxid	name
   ```

**Behavior**:
- **GTDB identifiers** (not NCBI) - no NCBIRanks needed
- Lineages reconstructed from report file indentation
- Tracks canonical ranks by depth
- Removes universal root (taxid=1) if present
- Strict validation: raises ValueError on format changes to catch version incompatibilities

**Edge cases**:
- Header presence/absence (optional since commit a17debb)
- Indentation must match rank hierarchy
- Non-canonical ranks validated to have canonical parents

### MetaMaps (`taxconverter/metamaps.py`)

```python
def parse_metamaps_krona(path: Path, lines: Iterable[str], ncbi: NCBIRanks) -> list[NCBIAnnotation]
```

**Input format**: Krona format, 3 columns
```
contig_name	taxID	name
```

**Behavior**:
- Skips comment lines (starting with `#`)
- taxID=0 means unclassified
- Uses NCBIRanks for lineage resolution

### MMseqs2 (`taxconverter/mmseqs.py`)

```python
def parse_mmseqs_tsv(path: Path) -> list[GenericAnnotation]
```

**Input format**: 9 tab-separated columns
- Column 9 contains semicolon-separated lineage (e.g., `d__Bacteria;p__Proteobacteria;...`)

**Behavior**:
- **GTDB identifiers** - no NCBIRanks needed
- Direct parsing from column 9
- Assumes lineage generated with `--tax-lineage` flag

## Output Formats

### Vamb Format (Default)

```
contigs	predictions
contig1	Bacteria;Proteobacteria;Gammaproteobacteria;...
contig2	unknown
```

- Header: `contigs\tpredictions`
- Lineage: semicolon-separated clade names
- Unclassified contigs get custom label (default: `unknown`, configurable with `--unassigned`)

### MMseqs Format (`--mmseqs-format` flag)

```
contig1	0	species	Escherichia	0	0	0	0	Bacteria;Proteobacteria;...;Escherichia
contig2	0	no rank	unknown	0	0	0	0	unknown
```

9 columns: `name\t0\trank\tlast_clade\t0\t0\t0\t0\tlineage`
- Rank determined by lineage length (0-6 = canonical, 7+ = subspecies, 0 = no rank)

## CLI Usage

```bash
# Centrifuge
taxconverter centrifuge -i input.tsv -o output.tsv

# Kraken2
taxconverter kraken2 -i input.tsv -o output.tsv

# Metabuli (requires TWO files)
taxconverter metabuli -c classifications.tsv -r report.tsv -o output.tsv

# MetaMaps
taxconverter metamaps -i input.krona -o output.tsv

# MMseqs2
taxconverter mmseqs2 -i input.tsv -o output.tsv

# Output to stdout
taxconverter centrifuge -i input.tsv

# MMseqs format output
taxconverter kraken2 -i input.tsv -o output.tsv --mmseqs-format

# Custom unassigned label
taxconverter centrifuge -i input.tsv -o output.tsv --unassigned "unclassified"
```

## Testing

### Test Data Locations

- **Small samples**: `test/data/centrifuge.tsv` (~20 lines)
- **Large samples**: `test/big_data/` (symlink to larger datasets)
  - `centrifuge_gi.tsv` (81k lines)
  - `kraken2_urog.tsv` (58k lines)
  - `Oral.metabuli_classifications.tsv` (202k lines)
  - `Oral.metabuli_report.tsv` (1.4k lines)
  - `oral_taxonomy_mmseq.tsv` (202k lines)
  - `METAMAPS_classification_results.EM.reads2Taxon.krona` (777 lines)

### Running Tests

```bash
pytest                                    # Run all tests
pytest tests/test_centrifuge.py -v       # Specific parser
pytest --cov=taxconverter                # With coverage
```

Tests configured in `pyproject.toml` under `[tool.pytest.ini_options]`.

## Development

### Installation

```bash
pip install -e .
```

Requires Python <3.12 (pandas 3.0 compatibility).

### Code Quality

Uses ruff for linting:
```bash
ruff check .
```

Config in `pyproject.toml` ignores E722 (bare except) and E501 (line too long).

### Version Management

Version defined in `taxconverter/__init__.py`:
```python
__version__ = (1, 1, 0)
```

Referenced dynamically in `pyproject.toml` via `tool.setuptools.dynamic`.

## Key Implementation Details

### Error Handling

- **Header validation**: Centrifuge requires exact header match
- **Field count**: All parsers validate field counts
- **NCBI ID lookup**: Raises `ValueError` if taxID not in NCBIRanks
- **Metabuli strictness**: Validates indentation, rank ordering, duplicate taxIDs
- **File existence**: CLI validates all input files exist before parsing
- **Output overwrite**: Raises `FileExistsError` if output exists

### Special Cases

- **taxID=0 or None**: Treated as unclassified (empty clades list)
- **Universal root (taxID=1)**: Metabuli removes if present
- **Non-canonical ranks**: Must have canonical parent (validated in NCBIRanks)
- **Empty lineages**: Get custom unassigned label in output
- **Semicolons**: Rejected in clade names (breaks output format)

### Performance Considerations

- NCBIRanks loads once at startup (~1-2s for compressed file)
- Uncompressed `.tsv` loads faster than `.tsv.gz` if available
- Session-scoped fixtures in tests avoid repeated loading
- Iterators used throughout for memory efficiency
