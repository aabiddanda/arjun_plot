"""Test ancillary plotting utilities."""
import pytest
import matplotlib.pyplot as plt
from arjun_plot import utils


@pytest.mark.parametrize(
    "test_in,test_seed,test_scale", [([-1, 0, 1], 42, 0.01)]
)  # noqa
def test_rand_jitter(test_in, test_seed, test_scale):
    """Testing random jitter."""
    jitter_res = utils.rand_jitter(test_in, scale=test_scale)
    assert test_in[0] != jitter_res[0]


def test_plot_yx():
    """Testing y=x plots."""
    _, ax = plt.subplots(1, 1)
    utils.plot_yx(ax)


def test_debox():
    """Testing of deboxing."""
    _, ax = plt.subplots(1, 1)
    utils.debox(ax)
    assert ax.spines["right"].get_visible() is False
    assert ax.spines["top"].get_visible() is False
    assert ax.spines["bottom"].get_visible()
    assert ax.spines["left"].get_visible()


def test_remove_border():
    """Testing removal of the border entirely."""
    _, ax = plt.subplots(1, 1)
    utils.remove_border(ax)
    for i in ax.spines:
        assert ax.spines[i].get_visible() is False


def test_remove_ticks():
    """Testing removal of all ticks."""
    _, ax = plt.subplots(1, 1)
    utils.remove_ticks(ax)
    assert ax.get_yticks().size == 0
    assert ax.get_xticks().size == 0


def test_label_multipanel():
    """Test for labeling of multipanel plots."""
    n = 3
    labels = ["A", "B", "C"]
    _, axs = plt.subplots(1, n)
    utils.label_multipanel(axs, labels)
    with pytest.raises(AssertionError):
        _, axs = plt.subplots(1, 2)
        utils.label_multipanel(axs, labels)
