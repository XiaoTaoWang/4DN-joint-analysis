import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import umap, joblib, sys, os, glob
from umap.plot import _datashade_points, _themes, _select_font_color
from collections import Counter
from palettable.tableau import Tableau_20
import matplotlib.gridspec as gridspec

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

    cmap = Tableau_20.get_mpl_colormap()
    counts = Counter(labels)
    new_labels = -np.ones_like(labels)
    i = 0
    for k in counts:
        if k == -1:
            continue
        if counts[k] < min_cluster:
            continue
        mask = labels == k
        new_labels[mask] = i
        i += 1
    
    color_code = {
        -1: '#FFFFFF',
    }

    max_i = new_labels.max()
    for i in range(max_i+1):
        c = cmap(i/max_i)
        color_code[i] = matplotlib.colors.rgb2hex(c)
    
    return new_labels, color_code

cell = sys.argv[1]
queue = glob.glob(os.path.join('{0}-Umap'.format(cell), '*pkl'))
community_fils = glob.glob(os.path.join('community.parameter-tuning.{0}'.format(cell), '*pkl'))
min_cluster = 500
for q in queue:
    coords = joblib.load(q)
    outfolder, outpre = os.path.split(q)
    outpre = outpre.replace('.pkl', '')
    Indicator = '{0}.completed'.format(outpre)
    lockFile = '{0}.lock'.format(outpre)
    if os.path.exists(Indicator):
        continue
    if os.path.exists(lockFile):
        continue

    lock = open(lockFile, 'w')
    lock.close()

    for f in community_fils:
        labels = np.r_[joblib.load(f)]
        new_labels, color_code = assign_colors(labels, min_cluster)
        k_, res_ = os.path.split(f)[1].replace('.pkl', '').split('_')[1:]
        out_fig = os.path.join(outfolder, '{0}.{1}.{2}.jpg'.format(outpre, k_, res_))
        ax = points_plot(
            coords,
            labels=new_labels,
            color_key=color_code,
            show_legend=True
        )
        plt.savefig(out_fig, dpi=1000)
        plt.close()
    
    completed = open(Indicator, 'w')
    completed.close()

    os.remove(lockFile)

