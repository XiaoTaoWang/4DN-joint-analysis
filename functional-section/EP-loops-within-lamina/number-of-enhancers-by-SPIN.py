import joblib, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pyensembl import EnsemblRelease

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def load_enhancers_by_gene(infil):

    hg38_db = EnsemblRelease(94, species='human')
    hg38_db.download()
    hg38_db.index()

    D = defaultdict(dict)
    with open(infil, 'r') as source:
        source.readline()
        for line in source:
            parse = line.rstrip().split('\t')
            if parse[0] in ['chrM', 'chrY']:
                continue
            if parse[6] == 'no':
                continue

            tc, ts, te, gid, tpm = parse[:5]
            tpm = float(tpm)
            gene = hg38_db.gene_by_id(gid)
            key = (tc, gene.start, gene.end, gene.strand, gid, tpm)
            enhancers = []
            tmp = parse[7]
            if tmp == '.':
                D[key] = enhancers
            else:
                for e in tmp.split(';'):
                    ec, es, ee = e.split(',')[1:]
                    es, ee = int(es), int(ee)
                    enhancers.append((ec, es, ee))
                D[key] = enhancers
    
    return D

def parse_SPIN_states(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e, name = line.rstrip().split()
            s, e = int(s), int(e)
            D[(c, s)] = name
    
    return D

cell = 'H1ESC'
enhancer_by_gene = load_enhancers_by_gene('{0}.distal-enhancers.tsv'.format(cell))
SPIN = parse_SPIN_states('{0}.SPIN.hg38.bed'.format(cell))
names = ['Speckle', 'Interior_Act1', 'Interior_Act2', 'Interior_Act3',
        'Interior_Repr1', 'Interior_Repr2', 'Near_Lm1', 'Near_Lm2', 'Lamina']

num_by_states = defaultdict(list)
for key in enhancer_by_gene:
    tmp_en = enhancer_by_gene[key]
    gc, gs, ge, strand = key[:4]
    if strand == '+':
        pos = gs//25000*25000
    else:
        pos = ge//25000*25000
    
    if not (gc, pos) in SPIN:
        continue

    g_state = SPIN[(gc, pos)]
    num_by_states[g_state].append(len(tmp_en))
    
pool = []
for key in names:
    pool.append(num_by_states[key])

colors = ['#812c43', '#bd4c58', '#e2745b',
          '#f2b179', '#f5d9a1', '#f1f4af',
          '#c8f0ba', '#92d2b9', '#6356a3']
fig = plt.figure(figsize=(1.8, 1.2))
ax = fig.add_subplot(111)
bp = ax.boxplot(pool, positions=[0, 1, 2, 3, 4, 5, 6, 7, 8], sym='', widths=0.6, patch_artist=True, notch=False,
                whiskerprops={'linestyle':'--', 'linewidth':0.7})
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(colors[i])
    patch.set_linewidth(0.7)
for median in bp['medians']:
    median.set(color='k', linewidth=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#ax.set_ylabel('RNA-Seq log2(TPM + 1)', fontsize=6)
ax.set_ylabel('Number of linked enhancers', fontsize=6)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
#ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.tick_params(axis='x', bottom=False, top=False, labelbottom=False, labeltop=False)
ax.set_axisbelow(True)
plt.savefig('{0}.number_of_enhancers.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()
