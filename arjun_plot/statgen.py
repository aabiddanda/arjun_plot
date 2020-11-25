"""Common plotting utilities I have used in statgen applications by Arjun Biddanda."""

import numpy as np


def qqplot_pval(ax, pvals, log_trans=True, **kwargs):
    """Plot a QQ-plot based on the quantiles of the p-values obtained."""
    assert np.all(pvals >= 0) and np.all(pvals <= 1)
    m = pvals.size
    # Calculate the expected quantiles of the p-value distribution
    exp_q = np.arange(1, m + 1) / (m + 1.0)
    true_q = np.sort(pvals)
    if log_trans:
        exp_q = -np.log10(exp_q)
        true_q = -np.log10(true_q)
    # Generate the plot and return the axis
    ax.scatter(exp_q, true_q, **kwargs)
    return ax
