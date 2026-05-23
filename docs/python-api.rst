.. _sec_python_api:

API Documentation
=================

Utility Functions
-----------------

General-purpose matplotlib helpers for axis styling, annotation, and data transforms.

.. automodule:: arjun_plot.utils
   :members:
   :undoc-members:
   :show-inheritance:


Admixture / STRUCTURE
---------------------

Bar-chart visualizations for ADMIXTURE/STRUCTURE ancestry proportion matrices.
Population sorting and cross-K colour matching are handled automatically.

.. automodule:: arjun_plot.admixture
   :members:
   :undoc-members:
   :show-inheritance:


Principal Components Analysis
------------------------------

The :class:`~arjun_plot.pca.PCA` class ingests SmartPCA or PLINK eigenvector output
and provides scatter plots, percent-variance axis labels, and population medoid overlays.

.. automodule:: arjun_plot.pca
   :members:
   :undoc-members:
   :show-inheritance:


Statistical Genetics / GWAS
----------------------------

Genome-wide and regional association visualizations, including Manhattan plots,
QQ-plots, locus-zoom plots with LD coloring, gene tracks, and single-variant
genotype/phenotype plots.

.. automodule:: arjun_plot.statgen
   :members:
   :undoc-members:
   :show-inheritance:


Karyograms
----------

Whole-genome ideogram (karyogram) figures with support for overlaying painted
ancestry tracts or other interval-based genomic features.

.. automodule:: arjun_plot.karyograms
   :members:
   :undoc-members:
   :show-inheritance:


Spatial / Geo-indexed
---------------------

Utilities for spatial and geographically indexed data visualizations.

.. automodule:: arjun_plot.spatial
   :members:
   :undoc-members:
   :show-inheritance:
