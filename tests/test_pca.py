"""Test PCA functions."""
import pytest
import matplotlib.pyplot as plt
from arjun_plot.pca import PCA


def test_read_smartpca():
    """Read a test smartpca output."""
    pca = PCA()
    pca.read_smartpca(evec_file="tests/test_data/test_smartpca.evec")


def test_read_plinkpca():
    """Read a test plink pca output."""
    pass


def test_add_poplabels():
    """Testing adding some population labels."""
    pass


def test_add_meta_data():
    """Testing adding in some new metadata."""
    pass


def test_check_data():
    """Testing that intermediate data checks make sense."""
    pass


def test_calc_percent_variance():
    """Testing that variance explained calculations make sense."""
    pass


def test_extract_medoid_pops_all():
    """Test extraction of medoids per-population grouping."""
    pass


def test_pca_axis_labels():
    """Test creation of PCA axis labels."""
    pass


def test_plot_pca():
    """Test plotting of PCA."""
    pass
