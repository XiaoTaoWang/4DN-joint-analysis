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


coords = joblib.load('HFFc6.euclidean-65.pkl')
new_labels = np.r_[joblib.load('consensus-clusters.new_labels.pkl')]
colors = [
    '#33A02C', '#1F78B4', '#FDBF6F',
    '#FB9A99', '#E31A1C', '#B2DF8A'
]
color_code = dict(zip(range(6), colors))
color_code[-1] = '#FFFFFF'
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
ax.set_xlim(170, 510)
ax.set_ylim(550, 210)
plt.savefig(out_fig, dpi=1000)
plt.close()