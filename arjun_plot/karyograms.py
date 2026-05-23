import warnings

import matplotlib as mpl
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch
from arjun_plot.utils import remove_border, remove_ticks

_DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)]

# GRCh38/hg38 primary-assembly chromosome lengths (base-pairs).
_hg38_chrom_sizes = pl.DataFrame(
    {
        "chrom": [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
        "size": [
            248_956_422, 242_193_529, 198_295_559, 190_214_555, 181_538_259,
            170_805_979, 159_345_973, 145_138_636, 138_394_717, 133_797_422,
            135_086_622, 133_275_309, 114_364_328, 107_043_718, 101_991_189,
             90_338_345,  83_257_441,  80_373_285,  58_617_616,  64_444_167,
             46_709_983,  50_818_468, 156_040_895,  57_227_415,
        ],
    }
)


def _resolve_color_map(categories, color_map, cmap):
    """Return a category → colour dict, sampling from *cmap* when *color_map* is None."""
    if color_map is None:
        cm = mpl.colormaps[cmap]
        n = len(categories)
        return {cat: cm(i / max(n, 1)) for i, cat in enumerate(categories)}
    missing = [c for c in categories if c not in color_map]
    assert not missing, f"color_map is missing entries for: {missing}"
    return color_map


def create_ideogram(chrom_df=None, chroms=None, scaling_factor=1e6, **kwargs):
    """Create a whole-genome ideogram with one subplot row per chromosome.

    Draws a scaled rectangle for each chromosome using the lengths provided in
    ``chrom_df``.  All spines and ticks are removed so that ancestry tracts or
    other features can be overlaid cleanly.

    :param polars.DataFrame | None chrom_df: DataFrame with columns ``chrom``
        (e.g. ``"chr1"``) and ``size`` (length in base-pairs).  Defaults to
        GRCh38/hg38 primary-assembly lengths when ``None``.
    :param list | None chroms: Ordered list of chromosome names to plot (default:
        the 22 autosomes ``["chr1", ..., "chr22"]``).  Pass e.g.
        ``["chr1", ..., "chr22", "chrX"]`` to include the X chromosome.
    :param float scaling_factor: Divisor applied to base-pair lengths before
        plotting (default 1e6 → lengths in Mb).
    :returns: ``(fig, axs, m_size)`` — the figure, list of axes (one per
        chromosome in the order given by *chroms*), and the maximum chromosome
        length in base-pairs.
    :rtype: tuple
    """
    assert scaling_factor > 0
    if chrom_df is None:
        chrom_df = _hg38_chrom_sizes
    assert "chrom" in chrom_df.columns
    assert "size" in chrom_df.columns

    if chroms is None:
        chroms = _DEFAULT_CHROMS

    absent = [c for c in chroms if c not in chrom_df["chrom"].to_list()]
    if absent:
        raise ValueError(f"chroms not found in chrom_df: {absent}")

    fig, axs = plt.subplots(len(chroms), 1, **kwargs)
    if len(chroms) == 1:
        axs = [axs]

    m_size = 0
    for i, chrom in enumerate(chroms):
        chrom_len = chrom_df.filter(pl.col("chrom") == chrom)["size"].to_numpy()[0]
        if chrom_len > m_size:
            m_size = chrom_len
        remove_border(axs[i])
        remove_ticks(axs[i])
        axs[i].add_patch(
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
        axs[i].plot([0, chrom_len / scaling_factor], [1, 1], color="none")
    return fig, axs, m_size


def draw_regions(axs, df=None, chroms=None, scaling_factor=1e6, categories=["conservative"], **kwargs):
    """Shade genomic regions (e.g. Neanderthal deserts) on existing ideogram axes.

    Filters ``df`` to rows whose ``type`` column matches any entry in ``categories``
    and shades the corresponding intervals with :func:`matplotlib.axes.Axes.axvspan`.
    Chromosomes not present in *chroms* are skipped with a warning.

    :param list axs: List of matplotlib axes returned by :func:`create_ideogram`.
    :param polars.DataFrame df: DataFrame with columns ``chrom``, ``start``, ``end``
        (base-pair coordinates), and ``type`` (region category label).
    :param list | None chroms: Ordered list of chromosome names corresponding to
        *axs* (default: the 22 autosomes).  Must match the *chroms* used when the
        ideogram was created.
    :param float scaling_factor: Same divisor used when creating the ideogram
        (default 1e6 → Mb).  Must match the value passed to :func:`create_ideogram`
        so shaded regions land at the correct positions.
    :param list categories: Category labels in ``df["type"]`` to include (default
        ``["conservative"]``).
    :param kwargs: Extra keyword arguments forwarded to
        :func:`matplotlib.axes.Axes.axvspan` (e.g. ``alpha``, ``color``).
    :returns: The updated ``axs`` list.
    :rtype: list
    """
    if chroms is None:
        chroms = _DEFAULT_CHROMS
    assert len(axs) == len(chroms), (
        f"axs has {len(axs)} entries but chroms has {len(chroms)}"
    )
    chrom_to_idx = {chrom: i for i, chrom in enumerate(chroms)}

    assert df is not None
    filt_deserts = df.filter(pl.col("type").is_in(categories))
    for chrom, start, end in zip(
        filt_deserts["chrom"].to_numpy(),
        filt_deserts["start"].to_numpy(),
        filt_deserts["end"].to_numpy(),
    ):
        if chrom not in chrom_to_idx:
            warnings.warn(f"{chrom} is not in the provided chroms list, skipping!")
            continue
        axs[chrom_to_idx[chrom]].axvspan(start / scaling_factor, end / scaling_factor, **kwargs)
    return axs


def plot_tracts(
    features_df, axs, chroms=None, scaling_factor=1e6, feat_idx=1, n_features=3, **kwargs
):
    """Overlay ancestry tracts or other genomic interval features onto ideogram axes.

    Each interval in ``features_df`` is drawn as a filled rectangle occupying a
    horizontal band within its chromosome's axis.  Multiple feature tracks can be
    stacked by varying ``feat_idx`` across calls (1-indexed, must be ≤ ``n_features``).
    The function also labels each chromosome axis on the left and returns a
    :class:`matplotlib.patches.Patch` legend handle built from the supplied kwargs.

    :param polars.DataFrame features_df: DataFrame with columns ``chrom``, ``start``,
        and ``end`` (positions in base-pairs).
    :param list axs: List of matplotlib axes returned by :func:`create_ideogram`.
    :param list | None chroms: Ordered list of chromosome names corresponding to
        *axs* (default: the 22 autosomes).  Must match the *chroms* used when the
        ideogram was created.
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

    if chroms is None:
        chroms = _DEFAULT_CHROMS
    assert len(axs) == len(chroms), (
        f"axs has {len(axs)} entries but chroms has {len(chroms)}"
    )

    for i, chrom in enumerate(chroms):
        axs[i].set_ylabel(
            chrom, rotation=0, va="center", ha="right", fontsize=8, labelpad=-12
        )
        segments = features_df.filter(pl.col("chrom") == chrom)
        for s, e in zip(segments["start"].to_numpy(), segments["end"].to_numpy()):
            axs[i].add_patch(
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


def plot_meta_karyogram(
    features_df,
    chrom_df=None,
    chroms=None,
    category_col="category",
    color_map=None,
    cmap="tab10",
    scaling_factor=1e6,
    alpha=0.6,
    **kwargs,
):
    """Plot a barcode-style meta-karyogram aggregating many short genomic intervals.

    Renders one thin vertical tick per interval across the requested chromosomes,
    coloured by a category column (e.g. archaic ancestry type).  When thousands of
    segments from many individuals are overlaid the tick density reveals the
    genome-wide distribution at a glance — dense regions indicate enrichment, gaps
    indicate deserts.  The visual style follows Figure 3A of Biddanda et al. (2026
    bioRxiv), where Neanderthal, Denisovan, and Ghost archaic segments are overlaid
    on a whole-genome ideogram.

    Each unique category occupies an equal horizontal band of the chromosome axis
    (``1 / n_categories`` of the total height), and each segment is drawn as a filled
    rectangle spanning its actual genomic coordinates rather than a tick at the midpoint.
    Rectangles are rendered via :class:`~matplotlib.collections.PatchCollection` so the
    function remains fast for large datasets.

    :param polars.DataFrame features_df: DataFrame with columns ``chrom``, ``start``,
        ``end`` (base-pair coordinates), and a category column (see ``category_col``).
    :param polars.DataFrame | None chrom_df: Chromosome size DataFrame passed to
        :func:`create_ideogram` (requires columns ``chrom`` and ``size``).  Defaults
        to GRCh38/hg38 primary-assembly lengths when ``None``.
    :param list | None chroms: Ordered list of chromosome names to plot (default:
        the 22 autosomes).  Pass e.g. ``["chr1", ..., "chr22", "chrX"]`` to include
        the X chromosome.
    :param str category_col: Column in ``features_df`` whose unique values define the
        colour groups (default ``"category"``).
    :param dict | None color_map: Mapping of category label → matplotlib colour string.
        Every unique value in ``features_df[category_col]`` must be present.  If
        ``None``, colours are sampled from ``cmap``.
    :param str cmap: Name of a matplotlib colormap used to auto-assign colours when
        ``color_map`` is ``None`` (default ``"tab10"``).  Any registered matplotlib
        colormap is accepted.
    :param float scaling_factor: Divisor applied to base-pair coordinates before
        plotting (default 1e6 → Mb).
    :param float alpha: Opacity of rectangles (default 0.6).
    :param kwargs: Additional keyword arguments forwarded to :func:`create_ideogram`
        (e.g. ``figsize``).
    :returns: ``(fig, axs, legend_handles)`` — the figure, list of axes (one per
        chromosome), and a list of :class:`~matplotlib.patches.Patch` legend handles
        (one per category, ready to pass to ``ax.legend()``).
    :rtype: tuple
    """
    for c in ["chrom", "start", "end", category_col]:
        assert c in features_df.columns, f"Missing required column: '{c}'"
    assert scaling_factor > 0

    if chroms is None:
        chroms = _DEFAULT_CHROMS

    categories = sorted(features_df[category_col].unique().to_list())
    color_map = _resolve_color_map(categories, color_map, cmap)

    fig, axs, m_size = create_ideogram(chrom_df=chrom_df, chroms=chroms, **kwargs)
    max_mb = m_size / scaling_factor

    for i, chrom in enumerate(chroms):
        ax = axs[i]
        ax.set_ylabel(
            chrom, rotation=0, va="center", ha="right", fontsize=8, labelpad=-12
        )
        ax.set_xlim(0, max_mb)

        chrom_feats = features_df.filter(pl.col("chrom") == chrom)
        if len(chrom_feats) == 0:
            continue

        n_cats = len(color_map)
        for i, (category, color) in enumerate(color_map.items()):
            cat_feats = chrom_feats.filter(pl.col(category_col) == category)
            if len(cat_feats) == 0:
                continue
            y_lo = i / n_cats
            h = 1 / n_cats
            starts = cat_feats["start"].to_numpy() / scaling_factor
            ends = cat_feats["end"].to_numpy() / scaling_factor
            patches = [
                plt.Rectangle((s, y_lo), e - s, h) for s, e in zip(starts, ends)
            ]
            pc = PatchCollection(patches, facecolor=color, edgecolor="none", alpha=alpha)
            ax.add_collection(pc)

    legend_handles = [
        Patch(facecolor=color, edgecolor="none", label=cat)
        for cat, color in color_map.items()
    ]
    return fig, axs, legend_handles
