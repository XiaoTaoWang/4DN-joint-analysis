import numpy as np
import matplotlib.pyplot as plt
import umap, hdbscan, matplotlib, joblib, os, sys
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

def flip_features(ori_loops, ori_fea):

    loops = []
    fea = []
    for i in range(len(ori_loops)):
        c1, s1, e1, c2, s2, e2, label = ori_loops[i]
        arr = ori_fea[i]
        anchor1 = arr[:11]
        anchor2 = arr[11:]
        if anchor1.max() < anchor2.max():
            anchor1, anchor2 = anchor2, anchor1

        loops.append((c1, s1, e1, c2, s2, e2, label))
        fea.append(np.r_[anchor1, anchor2])
    
    fea = np.r_[fea]

    return loops, fea

cell = sys.argv[1]
# step 1
# parameter test
for n_n in [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
    Indicator = '{0}.{1}.completed'.format(cell, n_n)
    lockFile = '{0}.{1}.lock'.format(cell, n_n)
    if os.path.exists(Indicator):
        continue

    if os.path.exists(lockFile):
        continue

    lock = open(lockFile, 'w')
    lock.close()

    loops, fea = joblib.load('{0}.umap-input.pkl'.format(cell))
    loops, fea = flip_features(loops, fea)

    mapper = umap.UMAP(n_neighbors=n_n, min_dist=0, n_components=2, metric='euclidean').fit(fea)
    joblib.dump(mapper.embedding_, '{0}.euclidean-{1}.pkl'.format(cell, n_n), compress=('xz', 3))
    for m_s in [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
        for m_c in [400, 500, 600]:
            labels = hdbscan.HDBSCAN(min_samples=m_s, min_cluster_size=m_c).fit_predict(mapper.embedding_)
            out_fig = '{0}.euclidean-{1}.hdbscan{2}_{3}.png'.format(cell, n_n, m_s, m_c)
            fig = plt.figure(figsize=(12, 12))
            ax = points_plot(mapper.embedding_, labels=labels, show_legend=True)
            plt.savefig(out_fig, dpi=1000)
            plt.close()
    
    completed = open(Indicator, 'w')
    completed.close()

    os.remove(lockFile)
