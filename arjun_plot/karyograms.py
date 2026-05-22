import warnings

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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


def draw_regions(axs, df=None, categories=["conservative"], **kwargs):
    """Shade genomic regions (e.g. Neanderthal deserts) on existing ideogram axes.

    Filters ``df`` to rows whose ``type`` column matches any entry in ``categories``
    and shades the corresponding intervals with :func:`matplotlib.axes.Axes.axvspan`.
    Chromosomes that cannot be parsed (e.g. ``chrX``) are silently skipped with a
    warning.

    :param list axs: List of 22 matplotlib axes returned by :func:`create_ideogram`.
    :param polars.DataFrame df: DataFrame with columns ``chrom``, ``start``, ``end``
        (base-pair coordinates), and ``type`` (region category label).
    :param list categories: Category labels in ``df["type"]`` to include (default
        ``["conservative"]``).
    :param kwargs: Extra keyword arguments forwarded to
        :func:`matplotlib.axes.Axes.axvspan` (e.g. ``alpha``, ``color``).
    :returns: The updated ``axs`` list.
    :rtype: list
    """
    assert len(axs) == 22
    assert df is not None
    filt_deserts = df.filter(pl.col("type").is_in(categories))
    for chrom, start, end in zip(
        filt_deserts["chrom"].to_numpy(),
        filt_deserts["start"].to_numpy(),
        filt_deserts["end"].to_numpy(),
    ):
        try:
            axid = int(chrom[3:]) - 1
            axs[axid].axvspan(start / 1e6, end / 1e6, **kwargs)
        except ValueError:
            warnings.warn(f"{chrom} is not currently parseable!")
    return axs


def plot_tracts(
    features_df, axs, scaling_factor=1e6, feat_idx=1, n_features=3, **kwargs
):
    """Overlay ancestry tracts or other genomic interval features onto ideogram axes.

    Each interval in ``features_df`` is drawn as a filled rectangle occupying a
    horizontal band within its chromosome's axis.  Multiple feature tracks can be
    stacked by varying ``feat_idx`` across calls (1-indexed, must be ≤ ``n_features``).
    The function also labels each chromosome axis on the left and returns a
    :class:`matplotlib.patches.Patch` legend handle built from the supplied kwargs.

    :param polars.DataFrame features_df: DataFrame with columns ``chrom``, ``start``,
        and ``end`` (positions in base-pairs).
    :param list axs: List of 22 matplotlib axes returned by :func:`create_ideogram`.
    :param float scaling_factor: Same divisor used when creating the ideogram (default
        1e6 → Mb coordinates).
    :param int feat_idx: 1-based index of this track within the stacked band (default
        1).
    :param int n_features: Total number of stacked tracks (determines band height,
        default 3).
    :param kwargs: Keyword arguments forwarded to :class:`matplotlib.patches.Rectangle`
        and used to build the legend handle.  Must include ``facecolor``, ``edgecolor``,
        and ``label``.
    :returns: ``(axs, leg_elem)`` — the updated axes list and a legend patch handle.
    :rtype: tuple
    """
    for c in ["chrom", "start", "end"]:
        assert c in features_df.columns
    assert n_features > 0
    assert feat_idx <= n_features
    assert feat_idx > 0
    assert scaling_factor > 0
    for i in range(1, 23):
        axs[i - 1].set_ylabel(
            f"chr{i}", rotation=0, va="center", ha="right", fontsize=8, labelpad=-12
        )
        segments = features_df.filter(pl.col("chrom") == f"chr{i}")
        for s, e in zip(segments["start"].to_numpy(), segments["end"].to_numpy()):
            axs[i - 1].add_patch(
                plt.Rectangle(
                    (s / scaling_factor, (feat_idx - 1) / n_features),
                    (e - s) / scaling_factor,
                    1 / n_features,
                    **kwargs,
                )
            )
    leg_elem = Patch(
        facecolor=kwargs["facecolor"],
        edgecolor=kwargs["edgecolor"],
        label=kwargs["label"],
    )
    return axs, leg_elem
