"""Tests for karyogram plotting utilities."""

import pytest
import polars as pl
import matplotlib.pyplot as plt
from arjun_plot.karyograms import create_ideogram, plot_tracts


@pytest.fixture
def chrom_df():
    return pl.DataFrame(
        {
            "chrom": [f"chr{i}" for i in range(1, 23)],
            "size": [200000000 + i * 1000000 for i in range(22)],
        }
    )


def test_create_ideogram(chrom_df):
    """Test basic ideogram creation returns expected structure."""
    fig, axs, m_size = create_ideogram(chrom_df=chrom_df)
    assert m_size > 0
    assert len(axs) == 22
    plt.close(fig)


def test_create_ideogram_no_chrom_df():
    """Test that omitting chrom_df raises ValueError."""
    with pytest.raises(ValueError):
        create_ideogram(chrom_df=None)


def test_create_ideogram_missing_columns():
    """Test that a DataFrame missing required columns raises AssertionError."""
    bad_df = pl.DataFrame({"chrom": ["chr1"], "length": [1000000]})
    with pytest.raises(AssertionError):
        create_ideogram(chrom_df=bad_df)


def test_plot_tracts(chrom_df):
    """Test plotting genomic tracts onto an existing ideogram."""
    _, axs, _ = create_ideogram(chrom_df=chrom_df)
    features_df = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [1000000, 2000000],
            "end": [5000000, 8000000],
        }
    )
    plot_tracts(features_df, axs)
    plt.close("all")
