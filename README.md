[![CI](https://github.com/aabiddanda/arjun_plot/actions/workflows/ci.yml/badge.svg)](https://github.com/aabiddanda/arjun_plot/actions/workflows/ci.yml) [![Documentation](https://github.com/aabiddanda/arjun_plot/actions/workflows/documentation.yaml/badge.svg)](https://github.com/aabiddanda/arjun_plot/actions/workflows/documentation.yaml) [![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/) [![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# arjun_plot: customized plotting routines in matplotlib

## Rationale

Over the course of doing research, I've had to generate **many** plots. Typically these have been using custom `matplotlib` routines. The package here has specific classes and modules which help with some common routines in plotting of genetic data.

## Modules

- `utils`: small utility functions for `matplotlib` axes — removing spines, adding jitter, swarm-plot coordinate transforms, and labeling multi-panel figures.
- `admixture`: routines for plotting ADMIXTURE/STRUCTURE ancestry bar plots. Handles sorting individuals within populations and matching population columns across multiple K values. Built on top of the [dystruct](https://github.com/tyjo/dystruct) package.
- `pca`: a `PCA` class for reading SmartPCA or PLINK eigenvector/eigenvalue output, adding population labels, and generating scatter plots with optional medoid overlays, percent-variance axis labels, and axis rotation.
- `statgen`: statistical genetics visualizations including QQ-plots, Manhattan plots, locus-zoom plots (with optional LD coloring), gene-region tracks (via UCSC API), and single-variant genotype/phenotype plots.
- `karyograms`: whole-genome ideogram (karyogram) figures using chromosome length data; supports overlaying painted tract intervals.
- `spatial`: placeholder for spatial and geo-indexed plotting utilities.

## Installation

```bash
git clone https://github.com/aabiddanda/arjun_plot.git
cd arjun_plot/
pip install .
```

Or directly from GitHub:

```bash
pip install git+https://github.com/aabiddanda/arjun_plot.git@main
```

## Stylesheets

In addition to these plotting routines, there are two `matplotlib` stylesheets useful when re-generating plots for presentations (particularly in LaTeX). They target fonts commonly used in presentations — such as Fira Sans — to keep plot typography consistent with slide decks.
