import matplotlib, bisect, joblib
import numpy as np 
import matplotlib.pyplot as plt
from collections import Counter

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def bisect_search(p, ref):

    cache = set()
    if not p[0] in ref:
        return cache
    
    List = ref[p[0]]
    idx = max(0, bisect.bisect(List, p[1:])-1)
    for q in List[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add((p[0],)+q)
    
    return cache

def parse_bed(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            if not c in D:
                D[c] = []
            D[c].append((s, e))
    
    for c in D:
        D[c].sort()
    
    return D

h1_enhancers = parse_bed('H1ESC.enhancer-elements.bed')
h1_promoters = parse_bed('H1ESC.promoter-elements.bed')
hff_enhancers = parse_bed('HFFc6.enhancer-elements.bed')
hff_promoters = parse_bed('HFFc6.promoter-elements.bed')
h1_anchor_by_gene = joblib.load('H1ESC.housekeeping.interacting-loci.pkl')
hff_anchor_by_gene = joblib.load('HFFc6.housekeeping.interacting-loci.pkl')

h1_collect = []
for ID in h1_anchor_by_gene:
    h1_anchors = h1_anchor_by_gene[ID]
    for a in h1_anchors:
        e_cache = bisect_search(a, h1_enhancers)
        p_cache = bisect_search(a, h1_promoters)
        if len(e_cache) and len(p_cache):
            h1_collect.append('E+P')
        elif len(e_cache) and (not len(p_cache)):
            h1_collect.append('E-only')
        elif (not len(e_cache)) and len(p_cache):
            h1_collect.append('P-only')
        else:
            h1_collect.append('None')
h1_counts = Counter(h1_collect)
h1_total = sum(h1_counts.values())

hff_collect = []
for ID in hff_anchor_by_gene:
    hff_anchors = hff_anchor_by_gene[ID]
    for a in hff_anchors:
        e_cache = bisect_search(a, hff_enhancers)
        p_cache = bisect_search(a, hff_promoters)
        if len(e_cache) and len(p_cache):
            hff_collect.append('E+P')
        elif len(e_cache) and (not len(p_cache)):
            hff_collect.append('E-only')
        elif (not len(e_cache)) and len(p_cache):
            hff_collect.append('P-only')
        else:
            hff_collect.append('None')
hff_counts = Counter(hff_collect)
hff_total = sum(hff_counts.values())

category_names = ['E+P', 'E-only', 'P-only', 'None']
h1_arr = [h1_counts[k] for k in category_names]
hff_arr = [hff_counts[k] for k in category_names]
data = np.array([h1_arr, hff_arr])
data_cum = data.cumsum(axis=1)
colors = ['#7FC97F', '#BEAED4', '#FDC086', '#CCCCCC']

fig = plt.figure(figsize=(0.6, 1.3))
ax = fig.add_subplot(111)
handles = []
legend_labels = []
labels = ['H1ESC', 'HFFc6']
for i, (name, color) in enumerate(zip(category_names, colors)):
    heights = data[:, i]
    starts = data_cum[:, i] - heights
    b = ax.bar(labels, heights, bottom=starts, width=0.77, color=color)
    handles.append(b)
    legend_labels.append(name)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)

ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)

ax.legend(handles[::-1], legend_labels[::-1], frameon=True, fontsize=5)
ax.set_ylabel('Number of loop anchors linked to housekeeping genes', fontsize=6)
plt.savefig('stratify.svg', dpi=300, bbox_inches='tight')
plt.close()