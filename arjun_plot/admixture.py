"""
Copyright (C) 2017-2018 Tyler Joseph <tjoseph@cs.columbia.edu>.

This file is part of Dystruct.

Modified by Arjun Biddanda <abiddanda@uchicago.edu> to be
a library of python functions to plot bars.
"""

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy
import scipy.cluster.vq
import scipy.stats


def match_one_Q(QK, Qprv):
    """Match the columns (populations) of Q to Qprv."""
    npops = QK.shape[1]
    distances = np.zeros((npops, npops))
    for id1, row1 in enumerate(QK.T):
        for id2, row2 in enumerate(Qprv.T):
            distances[id1, id2] = np.square(row1 - row2).sum()
    distances[:, npops - 1] = 1e10
    assignments = [i for i in range(npops)]
    for i in range(npops):
        label = np.unravel_index(np.nanargmin(distances), (npops, npops))
        distances[label[0]] = np.array([np.nan for i in range(npops)])
        distances[:, label[1]] = np.array([np.nan for i in range(npops)])
        assignments[label[1]] = label[0]
    return assignments


def subset_Q(Q, labels, order, subset=None):
    """Subset Q matrix using a subset of labels."""
    if subset is not None:
        sublabels = []
        subids = []
        for idx, label in enumerate(labels):
            if label in subset:
                sublabels.append(label)
                subids.append(idx)

        labels = sublabels
        order = subset
        Q = Q[subids]
    return (Q, labels, order)


def plot_k(ax, Q, lbls, order, colors, subset=None, spacing=2, bar_width=1, **kwargs):
    """Plot a single run of ADMIXTURE/STRUCTURE."""
    # Subset if we need to
    Q, labels, order = subset_Q(Q, lbls, order, subset)

    z = np.zeros(Q.shape[1])
    groups = [[z for i in range(spacing)] for pop in order]
    labels_grouped = [["" for i in range(spacing)] for pop in order]
    for label, mixt in zip(lbls, Q):
        groups[order.index(label)].append(mixt)
        labels_grouped[order.index(label)].append(label)

    for idx, group in enumerate(groups):
        group_z = group[:spacing]
        group_r = group[spacing:]
        if len(group_r) <= 1:
            continue

        # group individuals by major population contributor
        group_r = np.array(group_r)
        ids = np.argmax(group_r, axis=1)
        group_r = group_r[np.argsort(ids)]

        # sort within each population cluster
        ids = np.argmax(group_r, axis=1)
        sort_by = scipy.stats.mode(ids)[0][0]
        for idy in ids:
            rows = group_r[idy == ids].tolist()
            rows.sort(key=lambda row: row[sort_by])
            group_r[idy == ids] = np.array(rows)
        groups[idx] = group_z + group_r.tolist()

    labels_ordered = []
    for group in labels_grouped:
        start_idx = np.argwhere(np.array(group) != "")[0]
        start_idx = int(start_idx)
        length = len(group[start_idx:])
        mid = max(int(length / 2) - 1, 0)
        filt_group = []
        for i, g in enumerate(group):
            if i != mid + start_idx:
                filt_group.append("")
            else:
                filt_group.append(g)
        for label in filt_group:
            labels_ordered.append(label)
    Q = np.zeros((Q.shape[0] + spacing * len(groups), Q.shape[1]))
    idx = 0
    for g in groups:
        for m in g:
            Q[idx] = m
            idx += 1

    # Changing things to a pandas dataframe?
    ancestry_matrix = Q
    ancestry_matrix = pd.DataFrame(ancestry_matrix)
    ancestry_matrix.plot.bar(
        ax=ax,
        stacked=True,
        legend=False,
        width=bar_width,
        linewidth=bar_width / 10,
        color=colors,
    )

    # Setting plotting parameters
    ax.set_yticklabels([])
    ax.set_xticklabels(labels_ordered, **kwargs)
    ax.tick_params(axis=u"both", which=u"both", length=0)
    ax.patch.set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    return (ax, ancestry_matrix)
