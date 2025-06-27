import bisect, random, joblib, os, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LogNorm
from collections import defaultdict
from palettable.colorbrewer.diverging import RdYlBu_7_r

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)


cell = 'H1ESC'
# adjust the cluster labels
labels = joblib.load('consensus-clusters.new_labels.pkl')
new_labels = -np.ones_like(labels)
label_map = {
    -1: -1,
    0: 4,
    1: 5,
    2: 3,
    3: 0,
    4: 1,
    5: 2,
}
for i in label_map:
    new_labels[np.where(labels==i)[0]] = label_map[i]
labels = new_labels
joblib.dump(labels, 'consensus-clusters.new_labels.reassigned.pkl')

loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]
reorg_loops = {}
for i in sorted(set(labels)):
    if i == -1:
        continue
    idx = np.where(labels==i)[0]
    reorg_loops[i] = [loops[j][:6] for j in idx]

# adjust the cluster labels for HMM results
score_pool = joblib.load('{0}.hmm-intermediate.by_umap_class.pkl'.format(cell))
new_pool = {}
for i in score_pool:
    new_pool[label_map[i]] = score_pool[i]
score_pool = new_pool
joblib.dump(score_pool, '{0}.hmm-intermediate.by_umap_class.reassigned.pkl'.format(cell))

# update the fold-enrichment plot
state_briefs = ['AP', 'WP', 'PP', 'SE', 'WE', 'Insulator', 'TT', 'TE', 'WT', 'Repressed']
cluster_labels = list(map(str, sorted(reorg_loops)))
top_n = 3
collect_pair_labels = []
for i in sorted(reorg_loops):
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    farr = score_pool[i]
    farr[farr==0] = farr[farr>0].min()

    # plot top n features
    sort_table = []
    for ii in range(len(farr)):
        for jj in range(len(farr)):
            sort_table.append((farr[ii,jj], (ii,jj)))
    sort_table.sort(reverse=True)
    for ii in range(top_n):
        tmp_pairs = tuple([state_briefs[sort_table[ii][1][0]], state_briefs[sort_table[ii][1][1]]])
        tmp_pairs = tuple(sorted(tmp_pairs))
        tmp_label = '-'.join(tmp_pairs)
        if not tmp_label in collect_pair_labels:
            collect_pair_labels.append(tmp_label)

collect_pair_labels = [
    'PP-PP', 'PP-Repressed', 'PP-WP', 'Insulator-Insulator', 'Insulator-WP', 'Insulator-PP',
    'AP-TT', 'TT-TT', 'AP-SE', 'SE-WP', 'SE-SE', 'AP-AP', 'AP-WP', 'AP-TE', 'TE-TT',
]
fc_matrix = np.zeros((len(collect_pair_labels), len(reorg_loops)))
for i in sorted(reorg_loops):
    farr = score_pool[i]
    farr[farr==0] = farr[farr>0].min()
    sort_table = []
    for ii in range(len(farr)):
        for jj in range(len(farr)):
            tmp_pairs = tuple([state_briefs[ii], state_briefs[jj]])
            tmp_pairs = tuple(sorted(tmp_pairs))
            tmp_label = '-'.join(tmp_pairs)
            if tmp_label in collect_pair_labels:
                idx = collect_pair_labels.index(tmp_label)
                fc_matrix[idx, i] = max(farr[ii,jj], farr[jj,ii])

fig = plt.figure(figsize=(1.8, 2.1))
ax = fig.add_subplot(111)
fc_norm = fc_matrix / fc_matrix.max(axis=0)
data = pd.DataFrame(fc_norm, columns=cluster_labels, index=collect_pair_labels)
cg = sns.heatmap(data, cmap=RdYlBu_7_r.get_mpl_colormap(), annot=fc_matrix, ax=ax,
                 annot_kws={'fontsize':5}, fmt='.1f',
                 yticklabels=True, xticklabels=True)
plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=6)
plt.setp(cg.xaxis.get_majorticklabels(), rotation=0, fontsize=6)
#cg.xaxis.set_ticks_position('top')
#cg.yaxis.set_ticks_position('right')
cg.set_clip_on(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(0.5)
#ax.spines['left'].set_linewidth(0.5)
ax.xaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
ax.yaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
plt.savefig('{0}.clusters-features.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()