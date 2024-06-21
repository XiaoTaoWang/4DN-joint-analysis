import joblib, bisect, matplotlib
import numpy as np
from collections import defaultdict, Counter
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

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

def parse_SPIN_states(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e, name = line.rstrip().split()
            s, e = int(s), int(e)
            D[(c, s)] = name
    
    return D

cell = 'H1ESC'
platform = 'all'
elements = parse_bed('{0}.enhancer-elements.bed'.format(cell))
anchor_by_gene = joblib.load('{0}.anchor-by-gene.{1}.pkl'.format(cell, platform))
row_info = sorted(anchor_by_gene)
SPIN = parse_SPIN_states('{0}.SPIN.hg38.bed'.format(cell))

by_states = defaultdict(list)
total_lamina_genes = set()
lamina_genes_enhancers = set()
enhancer_state_by_gene = defaultdict(set)
for key in row_info:
    transcripts = anchor_by_gene[key]
    tmp_en = set()
    for tc, ts, te in transcripts:
        distal_anchors = transcripts[(tc, ts, te)]
        for a in distal_anchors:
            cache = bisect_search(a, elements)
            tmp_en.update(cache)
    
    gc, gs, ge, strand = key[:4]
    if strand == '+':
        pos = gs//25000*25000
    else:
        pos = ge//25000*25000
    
    if not (gc, pos) in SPIN:
        continue
    
    g_state = SPIN[(gc, pos)]
    
    if g_state == 'Lamina':
        total_lamina_genes.add(key)
    
    for ec, es, ee in tmp_en:
        pos = (es + ee) // 2
        pos = pos//25000*25000
        if (ec, pos) in SPIN:
            by_states[g_state].append(SPIN[(ec, pos)])
            if g_state == 'Lamina':
                lamina_genes_enhancers.add(key)
                enhancer_state_by_gene[key[4]].add(SPIN[(ec, pos)])

by_state_counts = {}
for key in by_states:
    by_state_counts[key] = Counter(by_states[key])

names = ['Speckle', 'Interior_Act1', 'Interior_Act2', 'Interior_Act3',
        'Interior_Repr1', 'Interior_Repr2', 'Near_Lm1', 'Near_Lm2', 'Lamina']
colors = ['#812c43', '#bd4c58', '#e2745b',
          '#f2b179', '#f5d9a1', '#f1f4af',
          '#c8f0ba', '#92d2b9', '#6356a3']
data = []
for key in names:
    stats = np.zeros(len(names))
    for i in range(len(names)):
        if names[i] in by_state_counts[key]:
            stats[i] = by_state_counts[key][names[i]]
    
    data.append(stats)

data = np.r_[data]
total_num = data.sum(axis=1)
data = data.T / data.sum(axis=1)
data = data.T
data_cum = data.cumsum(axis=1)

fig = plt.figure(figsize=(2.5, 1.2))
gs = gridspec.GridSpec(1, 6, wspace=0.35, top=0.9, bottom=0.1, left=0.1,
                       right=0.9)
ax = fig.add_subplot(gs[0:4])

handles = []
legend_labels = []
for i, (name, color) in enumerate(zip(names, colors)):
    heights = data[:, i]
    starts = data_cum[:, i] - heights
    b = ax.bar(list(range(len(names))), heights, bottom=starts, width=0.77, color=color)
    handles.append(b)
    legend_labels.append(name)

ax.set_xticks(list(range(len(names))))
ax.set_xticklabels(names, rotation=90, ha='center', va='top', fontsize=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_ylim(0, 1)

ax.set_ylabel('Fraction of linked enhancers in each SPIN state', fontsize=6)
ax.set_xlabel('SPIN state of genes', fontsize=6)

ax.legend(handles[::-1], legend_labels[::-1], ncol=1, bbox_to_anchor=(1, -0.05, 1, 1), loc='lower left',
          labelspacing=0.52, frameon=False, fontsize=5)

plt.savefig('{0}.enhancer-spin-states.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()

joblib.dump([total_lamina_genes, lamina_genes_enhancers], '{0}.lamina-genes.pkl'.format(cell))

    
