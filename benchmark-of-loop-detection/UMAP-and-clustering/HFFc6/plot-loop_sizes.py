import numpy as np
import matplotlib, joblib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker as mticker
from collections import Counter, defaultdict

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)


cell = 'HFFc6'
min_cluster = 5000
colors = [
    '#33A02C', '#1F78B4', '#FDBF6F',
    '#FB9A99', '#E31A1C', '#B2DF8A',
]

cluster_labels = joblib.load('consensus-clusters.new_labels.pkl')
size_by_cluster = defaultdict(list)
loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]

for i, loop in zip(cluster_labels, loops):
    if i == -1:
        continue
    c1, s1, e1, c2, s2, e2 = loop[:6]
    p1 = (s1 + e1) // 2
    p2 = (s2 + e2) // 2
    size_by_cluster[i+1].append(abs(p2 - p1))

cluster_labels = sorted(size_by_cluster)
sizes_list = [np.log10(np.r_[size_by_cluster[i]]) for i in cluster_labels]

fig = plt.figure(figsize=(2.1, 1.5))
ax = fig.add_subplot(111)
bp = sns.violinplot(data=sizes_list, ax=ax, palette=colors, inner='box', linewidth=0.5)
for violin in ax.collections[::2]:
    violin.set_alpha(0.8)
for line in ax.lines:
    line.set_linewidth(0.5)  # Reduce box line width
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Loop Sizes (bp)', fontsize=6)
ax.set_xticklabels(list(map(str, cluster_labels)), fontsize=5)
ax.set_xlabel('Cluster ID', fontsize=6)

ax.spines['bottom'].set_linewidth(0.8)
ax.spines['left'].set_linewidth(0.8)
ax.xaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
ax.set_axisbelow(True)
plt.savefig('dist-{0}.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()