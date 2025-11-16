import numpy as np
import polars as pl


def create_ideogram(chrom_df=None, scaling_factor=1e6, **kwargs):
	"""
    Setup axes + rectangles for karyograms
	"""
	assert scaling_factor > 0
    if chrom_df is None:
        raise ValueError('Chromosome lengths need to be specified')
    else:
        assert "chrom" in chrom_df.columns
        assert "size" in chrom_df.columns
        # NOTE: need to orient to specific chroms.
        fig, axs = plt.subplots(22, 1, **kwargs)
        m_size = 0
        for i in range(1,23):
            l = chrom_df.filter(pl.col('chrom') == f'chr{i}')['size'].to_numpy()[0]
            if l >= m_size:
                m_size = l
            remove_border(axs[i-1]);
            remove_ticks(axs[i-1]);
            axs[i-1].add_patch(plt.Rectangle((0, 0), l/scaling_factor, 1, ls="-", lw=1, ec="black", fc="none"))
            axs[i-1].plot([0,l/scaling_factor], [1,1], color='none')
        return fig, axs, m_size


def plot_tracts(features_df, axs, scaling_factor=1e6, n_features=3, **kwargs):
	"""
    Plotting method for tracts
	"""
	for c in ["chrom", "start", "end"]:
		assert c in features_df.columns
    assert scaling_factor > 0
    pass
