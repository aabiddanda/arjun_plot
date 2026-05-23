"""Unit testing for common statgen plots."""

import pytest
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import uniform
from arjun_plot.statgen import (
    qqplot_pval,
    manhattan_plot,
    locus_plot,
    locuszoom_plot,
    overlap_interval,
    plot_null_snps,
    gene_plot,
    plot_gene_region_worker,
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


def test_plot_null_snps_discontinuity():
    """Test plot_null_snps with positions containing a gap larger than the threshold."""
    np.random.seed(42)
    _, ax = plt.subplots(1, 1)
    pos = np.concatenate([np.linspace(0, 1e6, 100), np.linspace(7e6, 8e6, 100)])
    pvals = np.random.uniform(0.5, 5.0, size=pos.size)
    plot_null_snps(ax, pos=pos, pvals=pvals, threshold=5e6)


def test_locuszoom_with_lead_variant():
    """Test locuszoom plot with a lead variant and no LD matrix."""
    _, ax = plt.subplots(1, 1)
    np.random.seed(42)
    n = 20
    pos = np.linspace(0.15, 0.65, n)
    chroms = np.array(["chr1"] * n)
    pvals = uniform.rvs(size=n)
    variants = np.array([f"var{i}" for i in range(n)])
    locuszoom_plot(
        ax,
        chroms=chroms,
        variants=variants,
        pos=pos,
        pvals=pvals,
        chrom="chr1",
        position_min=0.1,
        position_max=0.7,
        lead_variant=variants[10],
    )


def test_locuszoom_with_ld_matrix():
    """Test locuszoom plot with a lead variant and an LD matrix."""
    _, ax = plt.subplots(1, 1)
    np.random.seed(42)
    n = 10
    pos = np.linspace(0.15, 0.65, n)
    chroms = np.array(["chr1"] * n)
    pvals = uniform.rvs(size=n)
    variants = np.array([f"var{i:02d}" for i in range(n)])
    ld_matrix = np.eye(n)
    locuszoom_plot(
        ax,
        chroms=chroms,
        variants=variants,
        pos=pos,
        pvals=pvals,
        chrom="chr1",
        position_min=0.1,
        position_max=0.7,
        lead_variant=variants[5],
        ld_variant_ids=variants.copy(),
        ld_matrix=ld_matrix,
    )


def test_plot_gene_region_invalid_build():
    """Test that an invalid genome build raises ValueError."""
    _, ax = plt.subplots(1, 1)
    with pytest.raises(ValueError):
        plot_gene_region_worker(ax, build="hg17")


def test_plot_gene_region_invalid_track():
    """Test that an invalid UCSC track raises ValueError."""
    _, ax = plt.subplots(1, 1)
    with pytest.raises(ValueError):
        plot_gene_region_worker(ax, track="RefSeq")


def test_locus_plot_violinplot():
    """Test locus_plot using a violinplot instead of a boxplot."""
    _, ax = plt.subplots(1, 1)
    np.random.seed(42)
    geno = np.random.binomial(2, 0.4, size=200)
    pheno = geno * 0.1 + np.random.normal(size=200)
    _, ns, _ = locus_plot(ax, geno, pheno, boxplot=False)
    assert np.sum(ns) == 200


def test_locus_plot_few_genotypes():
    """Test locus_plot warning when fewer than 3 genotype classes are observed."""
    _, ax = plt.subplots(1, 1)
    np.random.seed(42)
    geno = np.zeros(50, dtype=int)
    pheno = np.random.normal(size=50)
    with pytest.warns(UserWarning):
        locus_plot(ax, geno, pheno)
