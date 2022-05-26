import pytest
from arjun_plot import utils


@pytest.mark.parametrize("test_in,test_seed,test_scale", [([-1, 0, 1], 42, 0.01)])
def test_rand_jitter(test_in, test_seed, test_scale):
    """Testing random jitter."""
    jitter_res = utils.rand_jitter(test_in, scale=test_scale)
    assert test_in[0] != jitter_res[0]
