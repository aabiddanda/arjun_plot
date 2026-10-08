"""Unit testing for common statgen plots."""

import pytest
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import uniform
import gzip
import json
import shutil
import sys
import types
from arjun_plot.statgen import (
    qqplot_pval,
    manhattan_plot,
    locus_plot,
    locuszoom_plot,
    overlap_interval,
    plot_null_snps,
    gene_plot,
    plot_gene_region_worker,
    fetch_gene_annotations,
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


def test_rescale_axis(gene_file):
    """Test rescaling of the x-axis."""
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    ax = gene_plot(ax=ax, gene_file=gene_file)
    rescale_axis(ax)


@pytest.fixture
def gene_file():
    """Small UCSC genePred (ncbiRefSeq-style) table for offline tests."""
    return "tests/test_genes/test_ncbiRefSeq.txt"


def test_fetch_genes_local_txt(gene_file):
    """Test that a local genePred table is filtered to the chromosome & window."""
    records = fetch_gene_annotations(
        chrom="chr1", position_min=1000000, position_max=2000000, gene_file=gene_file
    )
    names = {g["name2"] for g in records}
    assert names == {"GENEA", "GENEB", "LOC1001"}
    assert len(records) == 4
    assert all(isinstance(g["txStart"], int) for g in records)


def test_fetch_genes_local_txt_gz(gene_file, tmp_path):
    """Test that gzipped genePred tables are read correctly."""
    gz_file = tmp_path / "test_ncbiRefSeq.txt.gz"
    with open(gene_file, "rb") as f_in, gzip.open(gz_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    records = fetch_gene_annotations(
        chrom="chr2", position_min=1000000, position_max=2000000, gene_file=gz_file
    )
    assert [g["name2"] for g in records] == ["GENEC"]


@pytest.mark.parametrize("as_api_response", [True, False])
def test_fetch_genes_local_json(gene_file, tmp_path, as_api_response):
    """Test reading a saved UCSC API response (or bare list of records)."""
    records = fetch_gene_annotations(
        chrom="chr1", position_min=0, position_max=int(1e7), gene_file=gene_file
    )
    data = {"ncbiRefSeq": records} if as_api_response else records
    json_file = tmp_path / "chr1.ncbiRefSeq.json"
    json_file.write_text(json.dumps(data))
    json_records = fetch_gene_annotations(
        chrom="chr1",
        position_min=1000000,
        position_max=2000000,
        gene_file=json_file,
    )
    assert {g["name2"] for g in json_records} == {"GENEA", "GENEB", "LOC1001"}


def test_fetch_genes_missing_file(tmp_path):
    """Test that a missing annotation file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        fetch_gene_annotations(gene_file=tmp_path / "missing.txt")


def test_fetch_genes_unsupported_format(tmp_path):
    """Test that an unrecognized file extension raises ValueError."""
    bad_file = tmp_path / "genes.gtf"
    bad_file.write_text("")
    with pytest.raises(ValueError):
        fetch_gene_annotations(gene_file=bad_file)


def test_fetch_genes_header_and_blank_lines(gene_file, tmp_path):
    """Test that header/comment lines and blank lines are skipped."""
    txt_file = tmp_path / "with_header.txt"
    with open(gene_file) as fp:
        body = fp.read()
    txt_file.write_text("#bin\tname\tchrom\n\n" + body + "\n")
    records = fetch_gene_annotations(
        chrom="chr1", position_min=1000000, position_max=2000000, gene_file=txt_file
    )
    assert {g["name2"] for g in records} == {"GENEA", "GENEB", "LOC1001"}


def test_fetch_genes_custom_genepred_cols(gene_file, tmp_path):
    """Test reading a genePred table without the leading bin column."""
    no_bin_file = tmp_path / "no_bin.tsv"
    with open(gene_file) as fp:
        lines = [line.split("\t", 1)[1] for line in fp]
    no_bin_file.write_text("".join(lines))
    cols = (
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames",
    )
    records = fetch_gene_annotations(
        chrom="chr1",
        position_min=1000000,
        position_max=2000000,
        gene_file=no_bin_file,
        genepred_cols=cols,
    )
    assert {g["name2"] for g in records} == {"GENEA", "GENEB", "LOC1001"}
    assert "bin" not in records[0]
    # The default columns would misalign this table and find nothing on chr1
    assert (
        fetch_gene_annotations(
            chrom="chr1",
            position_min=1000000,
            position_max=2000000,
            gene_file=no_bin_file,
        )
        == []
    )


def test_fetch_genes_window_boundaries(gene_file):
    """Test that genes touching the window edges are kept (inclusive overlap)."""
    # GENEB spans 1300000-1400000
    assert [
        g["name2"]
        for g in fetch_gene_annotations(
            chrom="chr1",
            position_min=1400000,
            position_max=1450000,
            gene_file=gene_file,
        )
    ] == ["GENEB"]
    assert [
        g["name2"]
        for g in fetch_gene_annotations(
            chrom="chr1",
            position_min=1250000,
            position_max=1300000,
            gene_file=gene_file,
        )
    ] == ["GENEB"]
    assert (
        fetch_gene_annotations(
            chrom="chr1",
            position_min=1400001,
            position_max=1450000,
            gene_file=gene_file,
        )
        == []
    )


def test_fetch_genes_local_json_gz(gene_file, tmp_path):
    """Test reading a gzipped UCSC API response."""
    records = fetch_gene_annotations(
        chrom="chr2", position_min=0, position_max=int(1e7), gene_file=gene_file
    )
    json_gz_file = tmp_path / "chr2.ncbiRefSeq.json.gz"
    with gzip.open(json_gz_file, "wt") as fp:
        json.dump({"ncbiRefSeq": records}, fp)
    json_records = fetch_gene_annotations(
        chrom="chr2",
        position_min=1000000,
        position_max=2000000,
        gene_file=json_gz_file,
    )
    assert [g["name2"] for g in json_records] == ["GENEC"]


class _FakeResponse:
    """Minimal stand-in for a requests.Response."""

    def __init__(self, payload, status_ok=True):
        self.payload = payload
        self.status_ok = status_ok

    def raise_for_status(self):
        if not self.status_ok:
            raise RuntimeError("HTTP error")

    def json(self):
        return self.payload


def _fake_requests(response, calls):
    """Build a fake requests module that records the URLs it is asked for."""

    def get(url, headers=None):
        calls.append(url)
        return response

    return types.SimpleNamespace(get=get)


def test_fetch_genes_remote_mocked(monkeypatch):
    """Test the remote UCSC query path without network access."""
    payload = {"knownGene": [{"name2": "GENEX", "txStart": 10, "txEnd": 20}]}
    calls = []
    monkeypatch.setitem(
        sys.modules, "requests", _fake_requests(_FakeResponse(payload), calls)
    )
    records = fetch_gene_annotations(
        build="hg19",
        track="knownGene",
        chrom="chr3",
        position_min=100,
        position_max=200,
    )
    assert records == payload["knownGene"]
    assert len(calls) == 1
    assert "genome=hg19;track=knownGene;chrom=chr3;start=100;end=200" in calls[0]


def test_fetch_genes_remote_http_error(monkeypatch):
    """Test that HTTP errors from the UCSC API are raised rather than parsed."""
    calls = []
    monkeypatch.setitem(
        sys.modules,
        "requests",
        _fake_requests(_FakeResponse({}, status_ok=False), calls),
    )
    with pytest.raises(RuntimeError):
        fetch_gene_annotations()


def test_gene_plot_local(gene_file):
    """Test plotting genes from a local file, keeping the longest transcript."""
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    ax = gene_plot(ax=ax, position_min=1e6, position_max=2e6, gene_file=gene_file)
    labels = sorted(t.get_text() for t in ax.texts)
    assert labels == ["GENEA→", "LOC1001→", "←GENEB"]


def test_gene_plot_local_name_filt(gene_file):
    """Test that name_filt removes matching genes when plotting from a local file."""
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    ax = gene_plot(
        ax=ax,
        position_min=1e6,
        position_max=2e6,
        gene_file=gene_file,
        name_filt=["^LOC"],
    )
    labels = sorted(t.get_text() for t in ax.texts)
    assert labels == ["GENEA→", "←GENEB"]


def test_gene_plot_local_without_requests(gene_file, monkeypatch):
    """Test that plotting from a local file does not require requests."""
    monkeypatch.setitem(sys.modules, "requests", None)
    _, ax = plt.subplots(1, 1, figsize=(6, 2))
    gene_plot(ax=ax, position_min=1e6, position_max=2e6, gene_file=gene_file)


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
