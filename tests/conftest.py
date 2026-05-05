"""Shared pytest fixtures for PAQR3 tests."""

from pathlib import Path

import pytest

# Root of the repository (one level above this file's parent).
_REPO_ROOT = Path(__file__).parent.parent
_DATA_DIR = _REPO_ROOT / "tests" / "files"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Path to the tests/data directory with small test files."""
    return _DATA_DIR


@pytest.fixture(scope="session")
def test_gtf(data_dir: Path) -> str:
    """Path to the small test GTF annotation file."""
    return str(data_dir / "new_pas_quant_test.gtf")


@pytest.fixture(scope="session")
def test_atlas(data_dir: Path) -> str:
    """Path to the small test PAS atlas BED file."""
    return str(data_dir / "atlas_test.bed")


@pytest.fixture(scope="session")
def test_bw_pos(data_dir: Path) -> str:
    """Path to the positive-strand test BigWig file."""
    return str(data_dir / "test_data.positive.bw")


@pytest.fixture(scope="session")
def test_bw_neg(data_dir: Path) -> str:
    """Path to the negative-strand test BigWig file."""
    return str(data_dir / "test_data.negative.bw")
