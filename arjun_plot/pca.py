"""Library of functions to plot PCA by Arjun Biddanda."""

import numpy as np
from matplotlib import cm
from adjustText import adjust_text
from scipy.optimize import curve_fit
import pandas as pd


class PCA:
    """Custom PCA plotting object."""

    def __init__(self):
        """Initialize the PCA object with empty fields."""
        self.indiv_labels = None
        self.pop_labels = None
        self.evecs = None
        self.evals = None
        self.meta_data = None

    def read_smartpca(self, evec_file, eval_file=None):
        """Read eigenvectors and eigenvalues in SmartPCA format.

        :param string evec_file: File with eigenvectors in smartPCA format.
        :param string eval_file: File with eigenvalues in smartpca format.

        """
        df = pd.read_csv(evec_file, sep=r"\s+")
        evecs = df.values[:, :-1]
        pop_labels = df.values[:, -1]
        indiv_labels = df.index.values
        if eval_file is None:
            evals = df.columns[1:].astype(np.float32)
        else:
            evals = np.loadtxt(eval_file)
        self.evals = evals
        self.evecs = evecs
        self.indiv_labels = indiv_labels
        self.pop_labels = pop_labels

    def read_plinkpca(self, evec_file, eval_file):
        """Read in PLINK-formatted eigenvectors and eigenvalues.

        :param string evec_file: File with eigenvectors in plink format.
        :param string eval_file: File with eigenvalues in plink format.
        """
        evals = np.loadtxt(eval_file)
        pcs = pd.read_csv(evec_file, sep=r"\s+", header=None, engine="python")
        evecs = pcs.values[:, 2:].astype(np.float32)
        indiv_labels = pcs.values[:, 0].astype(str)
        # Setting the underlying values after reading this
        self.evecs = evecs
        self.evals = evals
        self.indiv_labels = indiv_labels

    def add_poplabels(self, pop_dict, error=None):
        """Load population labels using a dictionary mapping indivID -> popID.

        :param dict pop_dict: dictionary mapping individual IDs to population IDs.
        """
        pop_labels = []
        for i in self.indiv_labels:
            try:
                pop_labels.append(pop_dict[i])
            except KeyError:
                if error is None:
                    raise KeyError(
                        f"Individual label {i} does not have a population label"
                    )
                else:
                    pop_labels.append(error)
        self.pop_labels = np.array(pop_labels).astype(str)

    def add_meta_data(self, meta_dict):
        """Add in metadata for each sample using a dict.

        The metadata contained by the dictionary can also be
        in a dictionary format with keys indicating language family, region, etc.

        :param dict meta_dict: dictionary mapping indivId -> {meta_dict}.

        """
        # Check that each individual has a key in the meta data
        keys = np.array([k for k in meta_dict.keys()])
        assert np.all(np.isin(self.indiv_labels, keys))
        self.meta_data = meta_dict

    def check_data(self):
        """Check dimensionality of data prior to running any plotting routines."""
        assert self.indiv_labels.size == self.pop_labels.size
        assert self.indiv_labels.size == self.evecs.shape[0]
        assert self.evecs.shape[1] >= 2
        if self.evals is not None:
            assert self.evals.size == self.evecs.shape[1]
        if self.meta_data is not None:
            keys = np.array([k for k in self.meta_data])
            assert np.all(np.isin(self.indiv_labels, keys))

    def calc_percent_var(self, pc1=1, pc2=2, extrapolate=False):
        """Calculate the proportion of variation explained by a PC.

        :param int pc1: first principal component.
        :param int pc2: second principal component.
        :param bool extrapolate: extrapolate if full-rank eigenvalues are not available.

        """
        pc1 = pc1 - 1
        pc2 = pc2 - 1
        assert pc1 != pc2
        assert (pc1 < self.evals.size) & (pc1 >= 0)
        assert (pc2 < self.evals.size) & (pc2 >= 0)
        a = {}
        if extrapolate:
            a.f = lambda x, a, b, c: a * np.exp(-b * x) + c
            popt, pcov = curve_fit(a.f, np.arange(1, self.evals.size + 1), self.evals)
            # Assuming that it is full-rank ...
            n = self.evecs.shape[0]
            sum_var_extrapolate = a.f(np.arange(1, n), *popt)
            sum_var_extrapolate[: self.evals.size] = self.evals
            sum_var = np.sum(sum_var_extrapolate)
        else:
            sum_var = np.sum(self.evals)
        # Calculating the proportion of variation explained
        prop_var = np.array(self.evals) / sum_var
        return (prop_var[pc1], prop_var[pc2])

    def pca_axis_labels(self, ax, pc1=1, pc2=2, extrapolate=False, **kwargs):
        """Plot PCA axis labels.

        :param matplotlib.pyplot.axis ax: axis to generate plot from.
        :param int pc1: first principal component  plotted.
        :param int pc2: second principal component plotted.
        :param bool extrapolate: extrapolate if full-rank eigenvalues are not available.

        """
        if self.evals is not None:
            vA, vB = self.calc_percent_var(pc1, pc2, extrapolate)
            ax.set_xlabel(r"PC%d (%0.2f%%)" % (pc1, vA * 100), **kwargs)
            ax.set_ylabel(r"PC%d (%0.2f%%)" % (pc2, vB * 100), **kwargs)
        else:
            ax.set_xlabel(r"PC%d" % pc1, **kwargs)
            ax.set_ylabel(r"PC%d" % pc2, **kwargs)
        return ax

    def extract_medoid_pops_all(self, max_pc=5, ctr_func=np.median):
        """Extract the medoid position across all PCs for a given population.

        :param int max_pc: maximum number of PCs to calculate population centroid.
        :param np.function ctr_func: centroid function for population values on a PC.

        """
        assert self.pop_labels is not None
        npcs = np.min([self.evals.size, max_pc])
        uniq_pops = np.unique(self.pop_labels)
        medoid_dict = {}
        for u in uniq_pops:
            idx = self.pop_labels == u
            ctrs = np.zeros(npcs)
            for i in range(npcs):
                cur_pc_pop = self.evecs[idx, i]
                ctrs[i] = ctr_func(cur_pc_pop)
            medoid_dict[u] = ctrs
        return medoid_dict

    def rotate_axes(self, ax, degree, pc1=1, pc2=2, **kwargs):
        """Rotating the axis for a PCA plot."""
        from matplotlib.transforms import Affine2D

        tr = Affine2D().rotate_deg(degree)
        ax = self.plot_pca(ax, pc1=pc1, pc2=pc2, transform=tr + ax.transData, **kwargs)
        return ax

    def plot_pca(
        self,
        ax,
        pc1=1,
        pc2=2,
        medoids=None,
        select_pops=None,
        text_lbls=True,
        colors=None,
        cmap="viridis",
        base_color="gray",
        legend=True,
        eps=0.05,
        **kwargs,
    ):
        """Full plotting method for PCA."""
        texts = []
        if select_pops is not None:
            assert len(select_pops) > 0
            if colors is not None:
                assert len(colors) >= len(select_pops)
            else:
                n = len(select_pops)
                right_colors = cm.get_cmap(cmap, n)
                colors = right_colors(np.linspace(0, 1, n))
            for i in range(len(select_pops)):
                cur_pop = select_pops[i]
                if medoids is not None:
                    cur_x, cur_y = (medoids[cur_pop][0], medoids[cur_pop][1])
                    # Plotting the scatter and the text alongside
                    ax.scatter(
                        cur_x, cur_y, color=colors[i], zorder=5, label=cur_pop, **kwargs
                    )
                    if not legend:
                        txt = ax.text(
                            cur_x + eps,
                            cur_y + eps,
                            cur_pop,
                            color=colors[i],
                            zorder=5,
                            snap=True,
                            ha="left",
                            va="bottom",
                            clip_on=True,
                        )
                        texts.append(txt)
                # The broader scatter-plot
                idx = self.pop_labels == cur_pop
                popx, popy = self.evecs[idx, pc1 - 1], self.evecs[idx, pc2 - 1]
                if medoids is None:
                    # Need to lower alpha here
                    ax.scatter(popx, popy, color=colors[i], zorder=5, **kwargs)
                else:
                    ax.scatter(popx, popy, color=colors[i], zorder=5, **kwargs)
            if len(texts) > 0:
                adjust_text(texts)

        # Lighter plotting here
        ax.scatter(
            self.evecs[:, pc1 - 1],
            self.evecs[:, pc2 - 1],
            color=base_color,
            zorder=-10,
            **kwargs,
        )
        return ax
