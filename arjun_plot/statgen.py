"""Common plotting utilities used in statgen applications."""

import numpy as np

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
    # TODO: should we have some option to thin p-values that are in the upper-levels?
    ax.scatter(exp_q, true_q, **kwargs)
    return ax


def manhattan_plot(
    ax,
    chroms,
    pos,
    pvals,
    chrom_def=HUMAN_CHROMS,
    thin=1,
    colors=["blue", "orange"],
    **kwargs,
):
    """Generate a Manhattan plot using some custom arguments.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        chroms (np.array): numpy array of chromosome values
        pos (np.array):  numpy array of positions (in basepairs).
        pvals (np.array): numpy array of p-values.
        chrom_def (list): list of all possible chromosomes.
        thin (int): argument to thin every i^th p-value. Set higher for faster plotting.
        colors (list): list of two colors to alternate between for chromosomes.

    Returns:
        ax (matplotlib.axis): axis containing the Manhattan plot.

    """
    assert chroms.size == pos.size
    assert pvals.size == pos.size
    assert len(colors) == 2
    assert thin >= 1

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
            ax.scatter(
                max_pos + cur_pos[::thin], cur_pvals[::thin], color=colors[i], **kwargs
            )
            xpos.append(max_pos + np.mean(cur_pos))
            max_pos += np.max(cur_pos)
            # increment the counter and mod to keep even / odd order
            i = (i + 1) % 2
    ax.set_xticks([])
    return (ax, xpos)
