import pytest
import matplotlib.pyplot as plt
from arjun_plot import utils


@pytest.mark.parametrize("test_in,test_seed,test_scale", [([-1, 0, 1], 42, 0.01)])
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
