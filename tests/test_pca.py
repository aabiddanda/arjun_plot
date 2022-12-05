"""Test PCA functions."""
import pytest
import numpy as np
from arjun_plot.pca import PCA


@pytest.fixture
def smartpca_evec():
    """Test smartpca eigenvector file."""
    return "tests/test_pca/test_smartpca.evec"


@pytest.fixture
def smartpca_eval():
    """Test smartpca eigenval file."""
    return "tests/test_pca/test_smartpca.eval"


@pytest.fixture
def base_pca_obj():
    """PCA object for testing manipulations."""
    pca = PCA()
    evecs = np.array([[0.1, 0.2, 0.3], [0.3, 0.05, 0.2], [-0.2, 0.1, -0.3]])
    evals = np.array([1.0, 0.5, 0.5])
    pca.evecs = evecs
    pca.evals = evals
    pca.indiv_labels = ["A", "B", "C"]
    pca.pop_labels = ["pop1", "pop2", "pop2"]
    return pca


@pytest.fixture
def good_pop_dict():
    """Good population label dictionary."""
    return {"A": "popX", "B": "popX", "C": "popY"}


@pytest.fixture
def incomplete_pop_dict():
    """Incomplete population dictionary."""
    return {"A": "popX", "B": "popX"}


@pytest.fixture
def incorrect_pop_dict():
    """Incorrect population dictionary."""
    return {"A": "popX", "B": "popX", "D": "popZ"}


def test_read_smartpca(smartpca_evec, smartpca_eval):
    """Read a test smartpca output."""
    pca = PCA()
    pca.read_smartpca(evec_file=smartpca_evec)
    pca.read_smartpca(evec_file=smartpca_evec, eval_file=smartpca_eval)
    with pytest.raises(FileNotFoundError):
        pca.read_smartpca(evec_file="")


def test_read_plinkpca():
    """Read a test plink pca output."""
    pass


def test_add_poplabels(
    base_pca_obj, good_pop_dict, incomplete_pop_dict, incorrect_pop_dict
):
    """Testing adding some population labels."""
    base_pca_obj.add_poplabels(good_pop_dict)
    with pytest.raises(KeyError):
        base_pca_obj.add_poplabels(incomplete_pop_dict)
        base_pca_obj.add_poplabels(incorrect_pop_dict)


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
