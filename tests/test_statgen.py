"""Unit testing for common statgen plots."""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import uniform
from arjun_plot.statgen import qqplot_pval, manhattan_plot, locus_plot


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


def test_locus_plot():
    """Test out the locus plot for analyses."""
    _, ax = plt.subplots(1, 1, figsize=(4, 4))
    np.random.seed(42)
    nsamples = 500
    geno = np.random.binomial(2, 0.05, size=nsamples)
    pheno = geno * 0.1 + np.random.normal(size=nsamples, scale=0.5)
    locus_plot(ax, geno, pheno)
