"""Unit testing for common statgen plots."""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import uniform
from arjun_plot.statgen import (
    qqplot_pval,
    manhattan_plot,
    locus_plot,
    locuszoom_plot,
    overlap_interval,
    gene_plot,
    rescale_axis,
)


def test_qqplot():
    """Test out the QQ-plot function."""
    _, ax = plt.subplots(1, 1)
    pvals = uniform.rvs(size=100)
    qqplot_pval(ax=ax, pvals=pvals, s=10)


def test_manhattan_plot():
    """Testing out the manhattan plot function."""
    _, ax = plt.subplots(1, 1, figsize=(8, 1))
    nsnps = 100
    chroms = []
    for c in range(1, 23):
        chroms.extend([f"chr{c}" for _ in range(nsnps)])
    chroms = np.array(chroms)
    pvals = uniform.rvs(size=chroms.size)
    pos = np.sort(uniform.rvs(size=pvals.size))
    manhattan_plot(ax, chroms=chroms, pos=pos, pvals=pvals)


def test_manhattan_plot_unsorted():
    """Testing out the manhattan plot function."""
    _, ax = plt.subplots(1, 1, figsize=(8, 1))
    nsnps = 1000
    chroms = []
    for c in range(1, 23):
        chroms.extend([f"chr{c}" for _ in range(nsnps)])
    chroms = np.array(chroms)
    pvals = uniform.rvs(size=chroms.size)
    pos = uniform.rvs(size=pvals.size)
    manhattan_plot(ax, chroms=chroms, pos=pos, pvals=pvals)


def test_locus_plot():
    """Test out the locus plot for analyses."""
    _, ax = plt.subplots(1, 1, figsize=(4, 4))
    np.random.seed(42)
    nsamples = 500
    geno = np.random.binomial(2, 0.1, size=nsamples)
    pheno = geno * 0.1 + np.random.normal(size=nsamples, scale=0.5)
    _, n, _ = locus_plot(ax, geno, pheno)
    assert np.sum(n) == nsamples


def test_locuszoom_plot():
    """Test out initial locuszoom plot."""
    _, ax = plt.subplots(1, 1, figsize=(8, 1))
    nsnps = 1000
    chroms = []
    for c in range(1, 23):
        chroms.extend([f"chr{c}" for _ in range(nsnps)])
    chroms = np.array(chroms)
    pvals = uniform.rvs(size=chroms.size)
    pos = uniform.rvs(size=pvals.size)
    variants = np.array([f"{c}:{p}" for (c, p) in zip(chroms, pvals)])
    locuszoom_plot(
        ax,
        chroms=chroms,
        variants=variants,
        pos=pos,
        pvals=pvals,
        chrom="chr1",
        position_min=1e-1,
        position_max=5e-1,
    )


def test_overlaps():
    """Test the overlap between two intervals."""
    assert overlap_interval((0, 5), (2, 6))
    assert not overlap_interval((0, 5), (5.1, 6))


def test_gene_plot():
    """Test the plotting of genes."""
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    gene_plot(ax=ax)


def test_rescale_axis():
    """Test rescaling of the x-axis."""
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    ax = gene_plot(ax=ax)
    rescale_axis(ax)
