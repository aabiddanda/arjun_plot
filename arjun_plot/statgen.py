"""Common plotting utilities used in statgen applications."""

import numpy as np
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
            cur_pos = cur_pos[::thin]
            cur_pvals = pvals[idx][::thin]
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
