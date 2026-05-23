Introduction
============

``arjun_plot`` is a personal library of ``matplotlib``-based plotting routines for statistical genetics. The primary goals are to make common plot types easy to produce inside Jupyter notebooks and to keep figures maximally configurable for publication and presentation use.

Modules
-------

**utils**
    Small axis-level helpers that reduce boilerplate in any matplotlib workflow:

    * ``debox`` / ``remove_border`` — strip top/right or all axis spines.
    * ``remove_ticks`` — hide tick marks for conceptual diagrams.
    * ``plot_yx`` — overlay a y = x reference line using the current axis limits.
    * ``label_multipanel`` — add panel labels (A, B, C, …) to a list of axes with configurable offsets.
    * ``rand_jitter`` — add proportional random noise to a 1-D array to reduce overplotting.
    * ``swarm`` — compute x-coordinates for a beeswarm / swarm plot without a heavy dependency.

**admixture**
    Bar-chart visualizations for ADMIXTURE/STRUCTURE ancestry proportion matrices (Q matrices).
    Individuals are sorted within populations by their dominant ancestry component, and population
    colours can be matched across multiple K values using ``match_one_Q``. Derives from the
    `dystruct <https://github.com/tyjo/dystruct>`_ package.

**pca**
    A ``PCA`` class that ingests SmartPCA (``.evec``/``.eval``) or PLINK-formatted eigenvector
    files and exposes methods for:

    * Adding population labels from a dictionary (``add_poplabels``).
    * Attaching arbitrary per-sample metadata (``add_meta_data``).
    * Computing the percent variance explained by any pair of PCs (``calc_percent_var``).
    * Generating publication-ready scatter plots with population medoid overlays,
      non-overlapping text labels, and optional axis rotation (``plot_pca``).

**statgen**
    Statistical genetics visualizations for genome-wide and regional analyses:

    * ``qqplot_pval`` — quantile–quantile plot for a set of p-values.
    * ``manhattan_plot`` — genome-wide Manhattan plot with alternating chromosome colours
      and convex-hull aggregation of null SNPs for compact vector output.
    * ``locuszoom_plot`` — regional association plot with optional LD-coloured scatter
      (R² computed against a lead variant) and grey shading for unmatched variants.
    * ``gene_plot`` / ``plot_gene_region_worker`` — gene/exon track fetched live from the
      UCSC Genome Browser API (hg19 or hg38).
    * ``locus_plot`` — boxplot or violinplot of phenotype distributions stratified by
      genotype at a single variant.
    * ``extract_ld_matrix`` — compute an empirical R² matrix from a tabix-indexed VCF
      over a specified genomic window.

**karyograms**
    Whole-genome ideogram (karyogram) figures. ``create_ideogram`` lays out 22 autosomal
    chromosomes as scaled rectangles given a data frame of chromosome lengths.
    ``plot_tracts`` overlays painted ancestry tracts or other interval-based features.

**spatial**
    Placeholder module for future spatial and geo-indexed data visualizations.

Notes
-----

This software is primarily a personal library. While it is reasonably well-tested, no
guarantee of bug-free operation is made. Bug reports and feature requests are welcome via
`GitHub Issues <https://github.com/aabiddanda/arjun_plot/issues>`_.
