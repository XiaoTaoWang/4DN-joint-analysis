import joblib, matplotlib, random
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

new_rc_params = {
    'text.usetex': False,
    'svg.fonttype': 'none'
}

matplotlib.rcParams.update(new_rc_params)

cell = 'HFFc6'
# load data
platforms = [
    'Hi-C', 'Micro-C', 'CTCF-ChIAPET', 'RNAPII-ChIAPET', 'H3K4me3-PLACSeq',
]
min_cluster = 5000
colors = [
    '#FB9A99', '#E31A1C', '#B2DF8A',
    '#FDBF6F', '#33A02C', '#1F78B4',
]
new_labels = joblib.load('consensus-clusters.new_labels.reassigned.pkl')
loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]
ratios = {}
for platform in platforms:
    pos_indices = []
    for i, loop in enumerate(loops):
        if platform in loop[-1]:
            pos_indices.append(i)
    tmp = []
    for i in pos_indices:
        if new_labels[i] > -1:
            tmp.append(new_labels[i]+1)
    counts = Counter(tmp)
    tmp = np.r_[[counts[i] for i in sorted(counts)]]
    ratios[platform] = tmp / tmp.sum()

fig = plt.figure(figsize=(1.1, 0.8))
ax = fig.add_subplot(111)
starts = np.zeros(len(platforms))
handles = []
for i, color in zip(np.arange(1,7), colors):
    tmp = np.r_[[ratios[platform][i-1] for platform in platforms]]
    b = ax.bar(platforms, tmp, bottom=starts, width=0.77, color=color)
    starts = starts + tmp
    handles.append(b)

legend_labels = list(map(str, np.arange(1,7)))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.8)
ax.spines['left'].set_linewidth(0.8)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Proporation', fontsize=6)
ax.set_xticklabels(platforms, ha='right', fontsize=5, rotation=45)
ax.legend(handles[::-1], legend_labels[::-1], ncol=1, bbox_to_anchor=(1.05, -0.02, 1, 1), loc='lower left',
            labelspacing=0.3, frameon=False, fontsize=5)

ax.xaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
ax.set_axisbelow(True)
plt.savefig('proporations-{0}.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()

        