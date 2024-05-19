# arjun_plot: customized plotting routines in matplotlib

## Rationale

Over the course of doing research, I've had to generate **many** plots. Typically these have been using custom `matplotlib` routines. The package here has specific classes and modules which help with some common routines in plotting of genetic data.

## Modules

- `utils`: contains mostly small utility functions for `matplotlib` axes
- `admixture`: contains routines for plotting ADMIXTURE bar-plots in (built off of the `dystruct` package)
- `pca`: some custom routines for plotting PCA in genetics applications
- `statgen`: custom plotting for statistical genetics applications (e.g. GWAS)

## Stylesheets

In addition to these plotting routines, I have generated two `matplotlib` stylesheets that I have found useful when re-generating plots for presentations (particularly in LaTeX). They try to match all of the fonts that I typically use for presentations - like Fira Sans
