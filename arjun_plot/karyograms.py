import polars as pl
import matplotlib.pyplot as plt
from arjun_plot.utils import remove_border, remove_ticks


def create_ideogram(chrom_df=None, scaling_factor=1e6, **kwargs):
    """Create a whole-genome ideogram with one subplot row per autosome.

    Draws a scaled rectangle for each of the 22 autosomes using the chromosome
    lengths provided in ``chrom_df``.  All spines and ticks are removed so that
    ancestry tracts or other features can be overlaid cleanly.

    :param polars.DataFrame chrom_df: DataFrame with columns ``chrom`` (e.g. ``"chr1"``) and ``size`` (length in base-pairs).
    :param float scaling_factor: Divisor applied to base-pair lengths before plotting (default 1e6 → lengths in Mb).
    :returns: ``(fig, axs, m_size)`` — the figure, array of axes (one per chromosome), and the maximum chromosome length in base-pairs.
    :rtype: tuple
    """
    assert scaling_factor > 0
    if chrom_df is None:
        raise ValueError("Chromosome lengths need to be specified")
    else:
        assert "chrom" in chrom_df.columns
        assert "size" in chrom_df.columns
        # NOTE: need to orient to specific chroms.
        fig, axs = plt.subplots(22, 1, **kwargs)
        m_size = 0
        for i in range(1, 23):
            chrom_len = chrom_df.filter(pl.col("chrom") == f"chr{i}")[
                "size"
            ].to_numpy()[0]
            if chrom_len >= m_size:
                m_size = chrom_len
                remove_border(axs[i - 1])
                remove_ticks(axs[i - 1])
                axs[i - 1].add_patch(
                    plt.Rectangle(
                        (0, 0),
                        chrom_len / scaling_factor,
                        1,
                        ls="-",
                        lw=1,
                        ec="black",
                        fc="none",
                    )
                )
                axs[i - 1].plot([0, chrom_len / scaling_factor], [1, 1], color="none")
    return fig, axs, m_size


def plot_tracts(features_df, axs, scaling_factor=1e6, n_features=3, **kwargs):
    """Overlay ancestry tracts or other genomic interval features onto ideogram axes.

    :param polars.DataFrame features_df: DataFrame with columns ``chrom``, ``start``, and ``end`` (positions in base-pairs).
    :param list axs: List of matplotlib axes returned by :func:`create_ideogram`.
    :param float scaling_factor: Same divisor used when creating the ideogram (default 1e6).
    :param int n_features: Maximum number of feature tracks to display per chromosome.
    """
    for c in ["chrom", "start", "end"]:
        assert c in features_df.columns
    assert scaling_factor > 0
    pass
