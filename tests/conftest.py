"""Shared pytest fixtures for taxconverter tests."""

import pytest
from pathlib import Path
from taxconverter.common import NCBIRanks


@pytest.fixture(scope="session")
def ncbi_ranks():
    """Load the NCBI taxonomy file (session-scoped for performance)."""
    # Uses the packaged data/clades.tsv file
    data_dir = Path(__file__).parent.parent / "data"
    clades_file = data_dir / "clades.tsv"

    # Try uncompressed first (faster), then compressed
    if clades_file.exists():
        return NCBIRanks.from_file(clades_file)
    else:
        clades_gz = data_dir / "clades.tsv.gz"
        if clades_gz.exists():
            return NCBIRanks.from_file(clades_gz)
        else:
            raise FileNotFoundError(
                f"Could not find NCBI lineage file at {clades_file} or {clades_gz}"
            )


@pytest.fixture
def fixtures_dir():
    """Path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def test_data_dir():
    """Path to existing test/data/ directory."""
    return Path(__file__).parent.parent / "test" / "data"
