import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import umap, joblib, sys, os, glob
from umap.plot import _datashade_points, _themes, _select_font_color
from collections import Counter
from palettable.colorbrewer.qualitative import Paired_12
import matplotlib.gridspec as gridspec

new_rc_params = {
    'text.usetex': False,
    'svg.fonttype': 'none'
}

matplotlib.rcParams.update(new_rc_params)

def points_plot(coords,
                labels=None,
                values=None,
                theme=None,
                cmap="Blues",
                color_key=None,
                color_key_cmap="Spectral",
                background="white",
                width=800,
                height=800,
                show_legend=True,
                subset_points=None,
                ax=None,
                alpha=None):
    
    if theme is not None:
        cmap = _themes[theme]["cmap"]
        color_key_cmap = _themes[theme]["color_key_cmap"]
        background = _themes[theme]["background"]

    if labels is not None and values is not None:
        raise ValueError(
            "Conflicting options; only one of labels or values should be set"
        )
    
    if alpha is not None:
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("Alpha must be between 0 and 1 inclusive")
    
    points = coords

    if subset_points is not None:
        if len(subset_points) != points.shape[0]:
            raise ValueError(
                "Size of subset points ({}) does not match number of input points ({})".format(
                    len(subset_points), points.shape[0]
                )
            )
        points = points[subset_points]

        if labels is not None:
            labels = labels[subset_points]
        if values is not None:
            values = values[subset_points]
    
    if points.shape[1] != 2:
        raise ValueError("Plotting is currently only implemented for 2D embeddings")
    
    font_color = _select_font_color(background)

    if ax is None:
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)
    
    if alpha is not None:
        alpha = alpha * 255
    else:
        alpha = 255

    ax = _datashade_points(
        points,
        ax,
        labels,
        values,
        cmap,
        color_key,
        color_key_cmap,
        background,
        width,
        height,
        show_legend,
        alpha,
    )

    ax.set(xticks=[], yticks=[])

    return ax

def assign_colors(labels, min_cluster):

    counts = Counter(labels)
    colors = [
        '#33A02C', '#1F78B4', '#FDBF6F',
        '#FB9A99', '#E31A1C', '#B2DF8A',
        '#A6CEE3'
    ]
    new_labels = -np.ones_like(labels)
    i = 0
    for k, count in counts.most_common():
        if k == -1:
            continue
        if count < min_cluster:
            continue
        mask = labels == k
        new_labels[mask] = i
        i += 1
    
    color_code = {
        -1: '#FFFFFF',
    }

    max_i = new_labels.max()
    for i in range(max_i+1):
        color_code[i] = colors[i]
    
    return new_labels, color_code

min_cluster = 5000
coords = joblib.load('H1ESC.euclidean-60.pkl')
labels = np.r_[joblib.load('consensus-clusters-5000_cutoff0.6_res0.5.pkl')]
new_labels, color_code = assign_colors(labels, min_cluster)
out_fig = 'umap.concensus.svg'
fig = plt.figure(figsize=(12, 12))
gs = gridspec.GridSpec(11, 11, wspace=0.05, top=0.95, bottom=0.05, left=0.05,
                    right=0.95)
ax = points_plot(
    coords,
    labels=new_labels,
    color_key=color_code,
    show_legend=True
)
ax.set_xlim(260, 540)
ax.set_ylim(560, 280)
plt.savefig(out_fig, dpi=1000)
plt.close()

joblib.dump(new_labels, 'consensus-clusters.new_labels.pkl')