"""Common plotting utilities used in statgen applications."""

import numpy as np
import warnings
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon


HUMAN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def qqplot_pval(ax, pvals, **kwargs):
    """Create a QQ-plot from a collection of pvalues.

    :param matplotlib.pyplot.axis ax: Input axis.
    :param np.array pvals: numpy array of p-values.

    """
    assert np.all(pvals >= 0) and np.all(pvals <= 1)
    m = pvals.size
    # Calculate the expected quantiles of the p-value distribution
    exp_q = np.arange(1, m + 1) / (m + 1.0)
    true_q = np.sort(pvals)
    # transform them to the negative log10 quantiles.
    exp_q = -np.log10(exp_q)
    true_q = -np.log10(true_q)
    # generate the plot as a scatter plot
    ax.scatter(exp_q, true_q, **kwargs)
    return ax


def plot_null_snps(ax, pos, pvals, q=0.80, threshold=5e6, **kwargs):
    """Make null-variants a convex hull polygon to save dots in vector graphics.

    :param matplotlib.pyplot.axis ax: Input axis.
    :param np.array pos:
    :param np.array pvals: negative log10 transformed p-values
    :param float q: quantile threshold for grouping effects
    :param float threshold: threshold for discontinuities in positioning
    """
    assert pos.size == pvals.size
    assert threshold > 1e6
    assert (q > 0) and (q < 1.0)
    # First determine if there is a discontinuity in SNPs due to the centromere ...
    discontinuities_idx = np.where(abs(np.diff(pos)) > threshold)[0] + 1
    if discontinuities_idx.size > 0:
        discontinuities_idx = np.insert(discontinuities_idx, 0, 0)
        discontinuities_idx = np.append(discontinuities_idx, pos.size - 1)
        for j, k in zip(discontinuities_idx[:-1], discontinuities_idx[1:]):
            cur_pos = pos[(pos >= pos[j]) & (pos < pos[k])]
            cur_pvals = pvals[(pos >= pos[j]) & (pos < pos[k])]

            null_pvals = cur_pvals[cur_pvals <= np.nanquantile(pvals, q)]
            null_pos = cur_pos[cur_pvals <= np.nanquantile(pvals, q)]
            assert null_pvals.size == null_pos.size
            if null_pvals.size > 0:
                hull_pts = np.vstack([null_pos, null_pvals]).T
                hull = ConvexHull(hull_pts)
                vertex_pts = np.vstack(
                    [hull_pts[hull.vertices, 0], hull_pts[hull.vertices, 1]]
                ).T
                ax.add_patch(Polygon(np.array(vertex_pts), **kwargs))
    else:
        null_pvals = pvals[pvals < np.nanquantile(pvals, q)]
        null_pos = pos[pvals < np.nanquantile(pvals, q)]
        assert null_pvals.size == null_pos.size
        if null_pvals.size > 0:
            hull_pts = np.vstack([null_pos, null_pvals]).T
            hull = ConvexHull(hull_pts)
            vertex_pts = np.vstack(
                [hull_pts[hull.vertices, 0], hull_pts[hull.vertices, 1]]
            ).T
            ax.add_patch(Polygon(np.array(vertex_pts), **kwargs))
    return ax, np.nanquantile(pvals, q)


def manhattan_plot(
    ax,
    chroms,
    pos,
    pvals,
    chrom_def=HUMAN_CHROMS,
    thin=1,
    colors=["blue", "orange"],
    q=0.80,
    threshold=10e6,
    padding=20e6,
    **kwargs,
):
    """Generate a Manhattan plot using some custom arguments.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        chroms (np.array): numpy array of chromosome values
        pos (np.array):  numpy array of positions (in basepairs).
        pvals (np.array): numpy array of p-values.
        chrom_def (list): list of all possible chromosomes.
        thin (int): thin every i^th p-value. Set higher for faster plotting.
        colors (list): list of two colors to alternate between for chromosomes.
        q (float): quantile below which to aggregate most p-values for plotting
        threshold (float): threshold for discontinuities
        padding (float): padding for between-chromosome distances
    Returns:
        ax (matplotlib.axis): axis containing the Manhattan plot.

    """
    assert chroms.size == pos.size
    assert pvals.size == pos.size
    assert len(colors) == 2
    assert thin >= 1
    assert padding > 0

    if np.all((pvals >= 0) & (pvals <= 1.0)):
        # Transform to the -log10 scale if necessary
        pvals = -np.log10(pvals)
    i = 0
    max_pos = 0
    xpos = []
    for x in chrom_def:
        idx = np.where(chroms == x)[0]
        if idx.size > 0:
            cur_pos = pos[idx] - np.min(pos[idx])
            cur_pvals = pvals[idx]
            # Sort the P-values by default
            pos_sort = np.argsort(cur_pos)
            cur_pvals = cur_pvals[pos_sort]
            cur_pos = cur_pos[pos_sort]
            cur_pos = cur_pos[::thin]
            cur_pvals = cur_pvals[::thin]
            plot_null_snps(
                ax,
                max_pos + cur_pos,
                cur_pvals,
                q=q,
                threshold=threshold,
                color=colors[i],
            )

            ax.scatter(
                max_pos + cur_pos[cur_pvals > np.nanquantile(cur_pvals, q)],
                cur_pvals[cur_pvals > np.nanquantile(cur_pvals, q)],
                color=colors[i],
                **kwargs,
            )
            xpos.append(max_pos + np.nanmax(cur_pos) / 2)
            # this padding value is kind of useful
            max_pos += np.nanmax(cur_pos) + padding
            # increment the counter and mod to keep even / odd order
            i = (i + 1) % 2
    # Set the xtick labels ...
    if len(chrom_def) > 0:
        ax.set_xticks([])
    return ax, xpos


def locuszoom_plot(
    ax, chroms, pos, pvals, chrom="chr1", position_min=1e6, position_max=2e6, **kwargs
):
    """Create a full locus-zoom plot for the GWAS summary statistics.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        chroms (np.array): numpy array of chromosome values
        pos (np.array):  numpy array of positions (in basepairs).
        pvals (np.array): numpy array of p-values.
        chrom (str): current chromosome being plotted.
        position_min (float): minimum position for a locus zoom
        position_max (float): maximum position for a locus zoom

    Returns:
        ax (matplotlib.axis): axis containing the locus-zoom plot
<<<<<<< HEAD
=======

>>>>>>> statgen
    """
    assert (position_min > 0) & (position_max > 0)
    assert position_max >= position_min
    assert chrom in chroms
    idx = (chroms == chrom) & (pos >= position_min) & (pos <= position_max)
    if np.all((pvals >= 0) & (pvals <= 1.0)):
        # Transform to the -log10 scale if necessary
        pvals = -np.log10(pvals)
    ax.scatter(pos[idx], pvals[idx], **kwargs)
    return ax


def overlap_interval(a, b):
    """Check if two intervals overlap at all.

    Args:
        a (tuple): A tuple for genomic range.
        b (tuple): A tuple for genomic range - can be None.

    Returns:
        overlap (bool): boolean indicator for overlapping.

    """
    if b is None:
        return False
    else:
        assert (len(a) == 2) & (len(b) == 2)
        assert (a[1] >= a[0]) & (b[1] >= b[0])
        return not ((a[1] < b[0]) or (b[1] < a[0]))


def overlap_labels(
    a, b, lbl_a, lbl_b, position_min, position_max, scaling_factor=0.015
):
    """Determine if two text labels sufficiently overlap.

    Args:
        a (tuple): A tuple for genomic range.
        b (tuple): A tuple for genomic range - can be None.
        lbl_a (string): Text label for interval a.
        lbl_b (string): Text label for interval b.
        position_min (float): minimum position for a locus zoom.
        position_max (float): maximum position for a locus zoom.
        scaling_factor (float): scaling factor to determine distance.

    Returns:
        overlap (bool): boolean indicator for overlapping text.

    """
    if (b is None) or (lbl_b is None):
        return False
    else:
        assert scaling_factor > 0.0
        position_range = position_max - position_min

        center_a = (a[0] + a[1]) / 2
        center_a_x = (center_a - position_min) / position_range

        a_width = scaling_factor * (1 + len(lbl_a))

        center_b = (b[0] + b[1]) / 2
        center_b_x = (center_b - position_min) / position_range
        b_width = scaling_factor * (1 + len(lbl_b))

        return overlap_interval(
            [center_a_x - a_width / 2, center_a_x + a_width / 2],
            [center_b_x - b_width / 2, center_b_x + b_width / 2],
        )


def plot_gene_region_worker(
    ax,
    build="hg38",
    chrom="chr1",
    position_min=1000000,
    position_max=1100000,
    yoff=0.1,
    scaling_factor=0.015,
    fontsize=6,
    name_filt=[],
):
    """Plot genes and exons.

    NOTE: much of the functionality is borrowed from https://github.com/krcurtis/locuszoom-plot/blob/master/locuszoom_plot/plot_gene_region.py

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        build (str): genome build and version.
        chrom (str): current chromosome being plotted.
        position_min (float): minimum position for a locus zoom.
        position_max (float): maximum position for a locus zoom.
        yoff (float): y offset for plotting gene names
        scaling_factor (float): scaling factor for text-labels (larger indicates less overlap).
        fontsize (float): fontsize of gene names.
        name_filt (list): list of regex elements to filter gene names by.

    Returns:
        ax (matplotlib.axis): axis containing the gene-region being plotted.

    """
    if build not in ["hg19", "hg38"]:
        raise ValueError(f"{build} is not a support genome build!")
    assert (position_min > 0) & (position_max > 0)
    assert position_max >= position_min
    assert yoff > 0
    assert scaling_factor > 0
    assert fontsize > 0
    import requests
    import re

    # Organize the gene-lists
    req = requests.get(
        f"https://api.genome.ucsc.edu/getData/track?genome={build};track=ncbiRefSeq;chrom={chrom};start={position_min};end={position_max}",  # noqa
        headers={"Content-Type": "application/json"},
    )
    results = req.json()["ncbiRefSeq"]
    names = []
    genes = []
    for g in results:
        # Want to avoid repeats
        if g["name2"] not in names:
            nmatch = 0
            for x in name_filt:
                nmatch += len(re.findall(x, g["name2"]))
            if nmatch == 0:
                genes.append(g)
                names.append(g["name2"])

    cur_x = None
    cur_lbl = None
    y_coord = 0
    for i, gene in enumerate(genes):
        center_x = (gene["txStart"] + gene["txEnd"]) / 2
        x = [float(gene["txStart"]), float(gene["txEnd"])]
        if not (
            overlap_interval(x, (position_min, position_min))
            or overlap_interval(x, (position_max, position_max))
        ):
            if overlap_interval(x, cur_x) or overlap_labels(
                x,
                cur_x,
                gene["name2"],
                cur_lbl,
                position_min=position_min,
                position_max=position_max,
                scaling_factor=scaling_factor,
            ):
                y_coord -= 0.5
            else:
                y_coord = 0.0
            cur_x = x
            y = [y_coord, y_coord]
            ax.plot(
                x, y, "-|", linewidth=1, color="black"
            )  ## force vertical line at beginning and end for short genes that might otherwis get dropped by the plot rendering engine
            exonStarts = [int(s) for s in gene["exonStarts"].split(",") if len(s) > 0]
            exonEnds = [int(s) for s in gene["exonEnds"].split(",") if len(s) > 0]
            for es, ee in zip(exonStarts, exonEnds):
                x = [es, ee]
                y = [y_coord, y_coord]
                ax.plot(x, y, "-|", linewidth=3, color="black")
            if "+" == gene["strand"]:
                ax.text(
                    center_x,
                    y_coord + yoff,
                    gene["name2"] + "→",
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fontsize,
                )
            elif "-" == gene["strand"]:
                ax.text(
                    center_x,
                    y_coord + yoff,
                    "←" + gene["name2"],
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fontsize,
                )
            cur_lbl = gene["name2"]
    return ax


def gene_plot(ax, chrom="chr1", position_min=1e6, position_max=2e6, **kwargs):
    """Plot the genes within a region."""
    ax = plot_gene_region_worker(
        ax=ax,
        chrom=chrom,
        position_min=int(np.round(position_min)),
        position_max=int(np.round(position_max)),
        **kwargs,
    )
    ax.set_yticks([])
    return ax


def rescale_axis(ax, scale=1e6):
    """Rescale the x-axis ticklabels to megabase."""
    assert scale > 1.0
    ax.set_xticklabels([f"{i}" for i in ax.get_xticks() / scale])
    return ax


def locus_plot(ax, genotypes, phenotypes, boxplot=True, **kwargs):
    """Plot the genotypes vs. phenotypes for a single-variant.

    Args:
        ax (matplotlib.axis): A matplotlib axis object to plot.
        genotypes (numpy.array): genotypes of each individual.
        phenotypes (numpy.array): phenotypes of each individual.
        boxplot (bool): display as a boxplot or violinplot.

    Returns:
        ax (matplotlib.axis): axis containing the variant-specific plot.
        ns (numpy.array): array of number of genotypes within each class.
        uniq_geno (numpy.array): array of the unique genotypes for phenotypes being plotted

    """
    assert genotypes.ndim == 1
    assert phenotypes.ndim == 1
    assert genotypes.size == phenotypes.size
    uniq_geno = np.sort(np.unique(genotypes))
    if uniq_geno.size < 3:
        warnings.warn("Less than 3 genotype classes observed!")
    ns = np.zeros(uniq_geno.size)
    for i, u in enumerate(uniq_geno):
        pheno = phenotypes[genotypes == u]
        ns[i] = pheno.size
        if boxplot:
            ax.boxplot(pheno, positions=[i], **kwargs)
        else:
            ax.violinplot(pheno, positions=[i], **kwargs)
    return ax, ns, uniq_geno

def extract_ld_matrix(vcf_fp, chrom='chr1', position_min=None, position_max=None):
    from cyvcf2 import VCF
    from tqdm import tqdm
    vcf =  VCF(vcf_fp, gts012=True, threads=12)
    ids = [] 
    gts = []
    for v in tqdm(vcf(f'{chrom}:{int(position_min)}-{int(position_max)}')):
        ids.append(v.ID)
        gts.append(v.gt_types.copy())
    assert len(ids) == len(gts)
    R2 = np.corrcoef(np.vstack(gts))**2
    ids = np.array(ids)
    assert R2.shape[0] == R2.shape[1]
    assert R2.shape[0] == ids.size
    return ids, R2
