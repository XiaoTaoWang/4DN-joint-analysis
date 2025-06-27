import joblib, matplotlib, random
import numpy as np
import matplotlib.pyplot as plt
from umap.plot import _datashade_points, _themes, _select_font_color

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

cell = 'HFFc6'
# load data
platforms = {'CTCF-ChIAPET', 'H3K4me3-PLACSeq', 'Hi-C',
             'Micro-C', 'RNAPII-ChIAPET'}
coords = joblib.load('{0}.euclidean-65.pkl'.format(cell))
loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]
for platform in platforms:
    pos_indices = []
    for i, loop in enumerate(loops):
        if platform in loop[-1]:
            pos_indices.append(i)

    neg_indices = []
    pos_indices_set = set(pos_indices)
    for i in range(len(coords)):
        if not i in pos_indices_set:
            neg_indices.append(i)
    
    if len(neg_indices) > len(pos_indices):
        neg_indices = random.sample(neg_indices, len(pos_indices))
    else:
        pos_indices = random.sample(pos_indices, len(neg_indices))
    indices = pos_indices + neg_indices
    sub = coords[indices]
    labels = -np.ones(len(sub), dtype=int)
    labels[:len(pos_indices)] = 1
    color_key = {1:'#E31A1C', -1:'#f1f1f1'}
    
    out_fig = '{0}.{1}.svg'.format(cell, platform)
    fig = plt.figure(figsize=(12, 12))
    ax = points_plot(
        sub,
        labels=labels,
        color_key=color_key,
        show_legend=False
    )
    ax.set_xlim(170, 510)
    ax.set_ylim(550, 210)
    plt.savefig(out_fig, dpi=1000)
    plt.close()