"""Common plotting utilities used in statgen applications."""

import numpy as np
import warnings
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon


HUMAN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def qqplot_pval(ax, pvals, **kwargs):
    """Create a QQ-plot from a collection of pvalues.

    :param matplotlib.pyplot.axis ax: Input axis.
    :param np.array pvals: numpy array of p-values.

    """
    assert np.all(pvals >= 0) and np.all(pvals <= 1)
    m = pvals.size
    # Calculate the expected quantiles of the p-value distribution
    exp_q = np.arange(1, m + 1) / (m + 1.0)
    true_q = np.sort(pvals)
    # transform them to the negative log10 quantiles.
    exp_q = -np.log10(exp_q)
    true_q = -np.log10(true_q)
    # generate the plot as a scatter plot
    ax.scatter(exp_q, true_q, **kwargs)
    return ax


def plot_null_snps(ax, pos, pvals, q=0.80, threshold=5e6, **kwargs):
    """Make null-variants a convex hull polygon to save dots in vector graphics.

    :param matplotlib.pyplot.axis ax: Input axis.
    :param np.array pos:
    :param np.array pvals: negative log10 transformed p-values
    :param float q: quantile threshold for grouping effects
    :param float threshold: threshold for discontinuities in variant positioning
    """
    assert pos.size == pvals.size
    assert threshold > 1e6
    assert (q > 0) and (q < 1.0)
    # First determine if there is a discontinuity in SNPs due to the centromere ...
    discontinuities_idx = np.where(abs(np.diff(pos)) > threshold)[0] + 1
    if discontinuities_idx.size > 0:
        discontinuities_idx = np.insert(discontinuities_idx, 0, 0)
        discontinuities_idx = np.append(discontinuities_idx, pos.size - 1)
        for j, k in zip(discontinuities_idx[:-1], discontinuities_idx[1:]):
            cur_pos = pos[(pos >= pos[j]) & (pos < pos[k])]
            cur_pvals = pvals[(pos >= pos[j]) & (pos < pos[k])]

            null_pvals = cur_pvals[cur_pvals <= np.nanquantile(pvals, q)]
            null_pos = cur_pos[cur_pvals <= np.nanquantile(pvals, q)]
            assert null_pvals.size == null_pos.size
            if null_pvals.size > 0:
                hull_pts = np.vstack([null_pos, null_pvals]).T
                hull = ConvexHull(hull_pts)
                vertex_pts = np.vstack(
                    [hull_pts[hull.vertices, 0], hull_pts[hull.vertices, 1]]
                ).T
                ax.add_patch(Polygon(np.array(vertex_pts), **kwargs))
    else:
        null_pvals = pvals[pvals < np.nanquantile(pvals, q)]
        null_pos = pos[pvals < np.nanquantile(pvals, q)]
        assert null_pvals.size == null_pos.size
        if null_pvals.size > 0:
            hull_pts = np.vstack([null_pos, null_pvals]).T
            hull = ConvexHull(hull_pts)
            vertex_pts = np.vstack(
                [hull_pts[hull.vertices, 0], hull_pts[hull.vertices, 1]]
            ).T
            ax.add_patch(Polygon(np.array(vertex_pts), **kwargs))
    return ax, np.nanquantile(pvals, q)


def manhattan_plot(
    ax,
    chroms,
    pos,
    pvals,
    chrom_def=HUMAN_CHROMS,
    thin=1,
    colors=["blue", "orange"],
    q=0.80,
    threshold=10e6,
    padding=20e6,
    **kwargs,
):
    """Generate a Manhattan plot using some custom arguments.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        chroms (np.array): numpy array of chromosome values
        pos (np.array):  numpy array of positions (in basepairs).
        pvals (np.array): numpy array of p-values.
        chrom_def (list): list of all possible chromosomes.
        thin (int): thin every i^th p-value. Set higher for faster plotting.
        colors (list): list of two colors to alternate between for chromosomes.
        q (float): quantile below which to aggregate most p-values for plotting
        threshold (float): threshold for discontinuities
        padding (float): padding for between-chromosome distances
    Returns:
        ax (matplotlib.axis): axis containing the Manhattan plot.

    """
    assert chroms.size == pos.size
    assert pvals.size == pos.size
    assert len(colors) == 2
    assert thin >= 1
    assert padding > 0

    if np.all((pvals >= 0) & (pvals <= 1.0)):
        # Transform to the -log10 scale if necessary
        pvals = -np.log10(pvals)
    i = 0
    max_pos = 0
    xpos = []
    for x in chrom_def:
        idx = np.where(chroms == x)[0]
        if idx.size > 0:
            cur_pos = pos[idx] - np.min(pos[idx])
            cur_pvals = pvals[idx]
            # Sort the P-values by default
            pos_sort = np.argsort(cur_pos)
            cur_pvals = cur_pvals[pos_sort]
            cur_pos = cur_pos[pos_sort]
            cur_pos = cur_pos[::thin]
            cur_pvals = cur_pvals[::thin]
            plot_null_snps(
                ax,
                max_pos + cur_pos,
                cur_pvals,
                q=q,
                threshold=threshold,
                color=colors[i],
            )

            ax.scatter(
                max_pos + cur_pos[cur_pvals > np.nanquantile(cur_pvals, q)],
                cur_pvals[cur_pvals > np.nanquantile(cur_pvals, q)],
                color=colors[i],
                **kwargs,
            )
            xpos.append(max_pos + np.nanmax(cur_pos) / 2)
            # this padding value is kind of useful
            max_pos += np.nanmax(cur_pos) + padding
            # increment the counter and mod to keep even / odd order
            i = (i + 1) % 2
    # Set the xtick labels ...
    ax.set_xticks([])
    return ax, xpos


def locuszoom_plot(
    ax, chroms, pos, pvals, chrom="chr1", position_min=1e6, position_max=2e6
):
    """Create a full locus-zoom plot."""
    raise NotImplementedError("Locuszoom plot is currently not supported!")


def plot_gene_region_worker(
    ax, build="b38", chrom="chr1", position_min=1e6, position_max=2e6
):
    """Helper function to plot genes and exons."""
    position_min = position_min / M
    position_max = position_max / M

    scale_gene_rows(gene_rows)

    text_yoffset = 0.4
    fontsize_magic = 5
    nrows = len(gene_rows)
    for i, row in enumerate(gene_rows):
        y_coord = -i
        for gene in row:
            center_x = (gene["txStart"] + gene["txEnd"]) / 2
            x = [float(gene["txStart"]), float(gene["txEnd"])]
            y = [y_coord, y_coord]
            ax.plot(
                x, y, "b-|", linewidth=1
            )  ## force vertical line at beginning and end for short genes that might otherwis get dropped by the plot rendering engine

            exonStarts = gene["exonStarts"]
            exonEnds = gene["exonEnds"]

            for es, ee in zip(exonStarts, exonEnds):
                x = [es, ee]
                y = [y_coord, y_coord]
                gene_axes.plot(x, y, "b-|", linewidth=5)

            if "+" == gene["strand"]:
                gene_axes.text(
                    center_x,
                    y_coord + text_yoffset,
                    gene["geneName"] + "→",
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fontsize_magic,
                )
            elif "-" == gene["strand"]:
                gene_axes.text(
                    center_x,
                    y_coord + text_yoffset,
                    "←" + gene["geneName"],
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fontsize_magic,
                )


def gene_plot(ax, chrom="chr1", start=1e6, end=2e6):
    """Plot the genes within a region."""
    raise NotImplementedError("Gene track plot is currently not implemented")


def locus_plot(ax, genotypes, phenotypes, boxplot=True, **kwargs):
    """Plot the genotypes vs. phenotypes for a single-variant.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        genotypes (numpy.array): genotypes of each individual
        phenotypes (numpy.array): phenotypes of each individual.
        boxplot (bool): display as a boxplot or violinplot
    Returns:
        ax (matplotlib.axis): axis containing the variant-specific plot.
        ns (numpy.array): array of number of genotypes within each class
        uniq_geno (numpy.array): array of the unique genotypes for phenotypes being plotted

    """
    assert genotypes.ndim == 1
    assert phenotypes.ndim == 1
    assert genotypes.size == phenotypes.size
    uniq_geno = np.sort(np.unique(genotypes))
    if uniq_geno.size < 3:
        warnings.warn("Less than 3 genotype classes observed!")
    ns = np.zeros(uniq_geno.size)
    for i, u in enumerate(uniq_geno):
        pheno = phenotypes[genotypes == u]
        ns[i] = pheno.size
        if boxplot:
            ax.boxplot(pheno, positions=[i], **kwargs)
        else:
            ax.violinplot(pheno, positions=[i], **kwargs)
    return ax, ns, uniq_geno
