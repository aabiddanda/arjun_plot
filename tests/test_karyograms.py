"""Tests for karyogram plotting utilities."""

import pytest
import polars as pl
import matplotlib.pyplot as plt
from arjun_plot.karyograms import create_ideogram, draw_regions, plot_tracts


@pytest.fixture
def chrom_df():
    return pl.DataFrame(
        {
            "chrom": [f"chr{i}" for i in range(1, 23)],
            "size": [200000000 + i * 1000000 for i in range(22)],
        }
    )


@pytest.fixture
def ideogram(chrom_df):
    fig, axs, m_size = create_ideogram(chrom_df=chrom_df)
    yield fig, axs, m_size
    plt.close(fig)


@pytest.fixture
def neanderthal_deserts():
    """Minimal Neanderthal desert regions spanning several autosomes."""
    return pl.DataFrame(
        {
            "chrom": ["chr3", "chr3", "chr7", "chr12", "chrX"],
            "start": [50_000_000, 120_000_000, 30_000_000, 80_000_000, 10_000_000],
            "end": [70_000_000, 140_000_000, 55_000_000, 95_000_000, 25_000_000],
            "type": [
                "conservative",
                "conservative",
                "vernot",
                "conservative",
                "conservative",
            ],
        }
    )


# ---------------------------------------------------------------------------
# create_ideogram
# ---------------------------------------------------------------------------


def test_create_ideogram(chrom_df):
    """Basic ideogram creation returns expected structure."""
    fig, axs, m_size = create_ideogram(chrom_df=chrom_df)
    assert m_size > 0
    assert len(axs) == 22
    plt.close(fig)


def test_create_ideogram_no_chrom_df():
    """Omitting chrom_df raises ValueError."""
    with pytest.raises(ValueError):
        create_ideogram(chrom_df=None)


def test_create_ideogram_missing_columns():
    """A DataFrame missing required columns raises AssertionError."""
    bad_df = pl.DataFrame({"chrom": ["chr1"], "length": [1000000]})
    with pytest.raises(AssertionError):
        create_ideogram(chrom_df=bad_df)


# ---------------------------------------------------------------------------
# draw_regions  (Neanderthal deserts test cases)
# ---------------------------------------------------------------------------


def test_draw_regions_conservative(ideogram, neanderthal_deserts):
    """Conservative deserts are shaded on the correct chromosome axes."""
    _, axs, _ = ideogram
    # axvspan adds Polygon patches; record baseline counts first
    baseline = [len(ax.patches) for ax in axs]
    result = draw_regions(
        axs,
        df=neanderthal_deserts,
        categories=["conservative"],
        alpha=0.3,
        color="steelblue",
    )
    assert result is axs
    # chr3 (index 2) has two conservative entries — patch count should grow by 2
    assert len(axs[2].patches) == baseline[2] + 2
    # chr12 (index 11) has one conservative entry
    assert len(axs[11].patches) == baseline[11] + 1


def test_draw_regions_vernot(ideogram, neanderthal_deserts):
    """Only vernot-category deserts are drawn when filtering by that type."""
    _, axs, _ = ideogram
    baseline = [len(ax.patches) for ax in axs]
    draw_regions(
        axs,
        df=neanderthal_deserts,
        categories=["vernot"],
        alpha=0.3,
        color="darkorange",
    )
    # chr7 (index 6) carries the vernot region — one new patch
    assert len(axs[6].patches) == baseline[6] + 1
    # chr3 should be untouched (only conservative entries there)
    assert len(axs[2].patches) == baseline[2]


def test_draw_regions_no_df(ideogram):
    """Passing df=None raises AssertionError."""
    _, axs, _ = ideogram
    with pytest.raises(AssertionError):
        draw_regions(axs, df=None)


def test_draw_regions_unparseable_chrom_warns(ideogram, neanderthal_deserts):
    """chrX entries produce a UserWarning and do not crash."""
    _, axs, _ = ideogram
    with pytest.warns(UserWarning):
        draw_regions(
            axs,
            df=neanderthal_deserts,
            categories=["conservative"],
            alpha=0.3,
            color="steelblue",
        )


# ---------------------------------------------------------------------------
# plot_tracts
# ---------------------------------------------------------------------------


def test_plot_tracts(ideogram):
    """Genomic tracts are drawn onto the ideogram without error."""
    _, axs, _ = ideogram
    features_df = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [1_000_000, 2_000_000],
            "end": [5_000_000, 8_000_000],
        }
    )
    axs, leg = plot_tracts(
        features_df,
        axs,
        facecolor="salmon",
        edgecolor="none",
        label="Neanderthal",
    )
    assert leg.get_label() == "Neanderthal"


def test_plot_tracts_stacked(ideogram):
    """feat_idx and n_features control the vertical position of tracks."""
    _, axs, _ = ideogram
    features_df = pl.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [10_000_000],
            "end": [20_000_000],
        }
    )
    # Two stacked tracks on the same ideogram should not raise
    plot_tracts(
        features_df,
        axs,
        feat_idx=1,
        n_features=2,
        facecolor="steelblue",
        edgecolor="none",
        label="track1",
    )
    plot_tracts(
        features_df,
        axs,
        feat_idx=2,
        n_features=2,
        facecolor="salmon",
        edgecolor="none",
        label="track2",
    )


def test_plot_tracts_missing_columns(ideogram):
    """A features DataFrame missing required columns raises AssertionError."""
    _, axs, _ = ideogram
    bad_df = pl.DataFrame({"chrom": ["chr1"], "begin": [0], "end": [1000]})
    with pytest.raises(AssertionError):
        plot_tracts(bad_df, axs, facecolor="blue", edgecolor="none", label="x")
