import numpy as np
import matplotlib.pyplot as plt
import joblib, matplotlib
from palettable.colorbrewer.sequential import Oranges_3, Purples_3
from adjustText import adjust_text

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

tf_pool = joblib.load('H1-union-loops.TF-enrich.pkl')
for key in tf_pool:
    Map = tf_pool[key]
    TFs = list(Map)
    ratio = [Map[tf][0] for tf in TFs]
    enrich = [Map[tf][1] for tf in TFs]
    enrich = np.r_[enrich]
    y = enrich
    x = np.r_[ratio]
    ratio_x = (x - x.mean()) / x.std()
    ratio_y = (y - y.mean()) / y.std()
    size = ratio_x + ratio_y

    TF_show = ['EZH2', 'POLR2A', 'CHD1', 'KDM4A', 'PHF8', 'TAF1', 'RAD21', 'CTCF']
    
    TF_show = set(TF_show)
    print('high ratio:')
    high_ratio = np.argsort(x)[::-1][:2]
    for i in high_ratio:
        print(TFs[i])
        TF_show.add(TFs[i])
    print('high enrich:')
    high_enrich = np.argsort(y)[::-1][:2]
    for i in high_enrich:
        print(TFs[i])
        TF_show.add(TFs[i])
    
    fig = plt.figure(figsize=(1, 0.9))
    ax = fig.add_subplot(111)
    ax.scatter(x[:len(TFs)], y[:len(TFs)], c='#FDAE6B', ec='none', s=9, alpha=0.9)
    
    texts = []
    for i in range(x.size):
        if TFs[i] in TF_show:
            texts.append(ax.text(x[i], y[i], TFs[i], va='center', ha='center', fontsize=5))
    adjust_text(texts, arrowprops=dict(arrowstyle='->, head_length=1.5, head_width=1.5', color='k', lw=0.5))

    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, 1)
    
    ax.set_xlabel('Fraction of loop anchors', fontsize=6)
    ax.set_ylabel('Fold Enrichment', fontsize=6)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.7)
    ax.spines['left'].set_linewidth(0.7)
    ax.xaxis.set_tick_params(width=0.7, pad=2, labelsize=6)
    ax.yaxis.set_tick_params(width=0.7, pad=2, labelsize=6)

    plt.savefig('cluster{0}.enrich-plot-loops.svg'.format(key), dpi = 300, bbox_inches='tight')
    plt.close()