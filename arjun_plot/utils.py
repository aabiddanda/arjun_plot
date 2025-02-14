"""Utility plotting functions in matplotlib."""
import numpy as np


def plot_yx(ax, **kwargs):
    """
    Plot the y=x line within a matplotlib axis.

    :param matplotlib.pyplot.axis ax: Input axis.
    """
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, **kwargs)


def debox(ax):
    """
    Remove the top and right spines of a plot.

    :param matplotlib.pyplot.axis ax: Input axis.
    """
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def remove_border(ax):
    """
    Remove border spines completely from a matplotlib axis.

    :param matplotlib.pyplot.axis ax: Input axis.
    """
    for i in ["top", "bottom", "left", "right"]:
        ax.spines[i].set_visible(False)


def remove_ticks(ax):
    """
    Remove the x and y ticks for conceptual plots.

    :param matplotlib.pyplot.axis ax: Input axis.
    """
    ax.set_yticks([])
    ax.set_xticks([])


def label_multipanel(axs, labels, xoff=-0.05, yoff=1.14, **kwargs):
    """
    Labeling multiple axes with text labels.

    :param list axs: Input axis.
    :param list labels: Input labels for each axis.
    :param float xoff: x-axis offset for axis label
    :param float yoff: y-axis offset for axis label

    """
    assert len(axs) == len(labels)
    for i, lbl in enumerate(labels):
        axs[i].text(xoff, yoff, lbl, transform=axs[i].transAxes, **kwargs)


def rand_jitter(arr, scale=0.01, seed=None):
    """
    Randomly jitter an array to avoid overlapping points.

    :param numpy.array arr: numpy array for applying jitter.
    :param float scale: scale for warping based on range in data.
    :param int seed: random seed.

    """
    assert scale > 0.0
    if seed is not None:
        assert type(seed) is int
        assert seed > 0
        np.random.seed(seed)
    stdev = scale * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev


def swarm(arr, nbins=None, width=1.0):
    """
    Swarm plot coordinate transform for ``arr``.
    Based on answer here: https://stackoverflow.com/a/76405274
    """
    assert width > 0
    y = np.asarray(arr)
    assert y.ndim == 1
    assert y.size > 0
    if nbins is None:
        nbins = np.ceil(y.size / 6).astype(int)

    # Get upper bounds of bins
    x = np.zeros(y.size)

    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()

    # Divide indices into bins
    ibs = []
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        i = np.nonzero((y > ymin) * (y <= ymax))[0]
        ibs.append(i)

    # Assign x indices
    dx = width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(yy)]
            a = i[j::2]
            b = i[(j + 1) :: 2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx
    return x
