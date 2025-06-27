import matplotlib
import numpy as np
from collections import defaultdict
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

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

def count_enhancers_within_lamina(cell):

    enhancer_by_gene = load_enhancers_by_gene('{0}.distal-enhancers.tsv'.format(cell))
    SPIN = parse_SPIN_states('{0}.SPIN.hg38.bed'.format(cell))
    active_gene_enhancers = []
    inactive_gene_enhancers = []
    for key in enhancer_by_gene:
        tmp_en = enhancer_by_gene[key]
        gc, gs, ge, strand, gid, tpm = key
        if strand == '+':
            pos = gs//25000*25000
        else:
            pos = ge//25000*25000
        
        if not (gc, pos) in SPIN:
            continue

        g_state = SPIN[(gc, pos)]
        if g_state == 'Lamina':
            if tpm > 1:
                active_gene_enhancers.append(len(tmp_en))
            else:
                inactive_gene_enhancers.append(len(tmp_en))
    
    return active_gene_enhancers, inactive_gene_enhancers

h1_active, h1_inactive = count_enhancers_within_lamina('H1ESC')
hff_active, hff_inactive = count_enhancers_within_lamina('HFFc6')

colors = ['#7FC97F', '#999999', '#BEAED4', '#999999']
fig = plt.figure(figsize=(1.4, 1.2))
ax = fig.add_subplot(111)
bp = ax.boxplot([h1_active, h1_inactive, hff_active, hff_inactive],
                positions=[0, 1, 2, 3], sym='', widths=0.6, patch_artist=True,
                notch=False, whiskerprops={'linestyle':'--', 'linewidth':0.7})
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(colors[i])
    patch.set_linewidth(0.7)
for median in bp['medians']:
    median.set(color='k', linewidth=1)

# pair 1
maxv = max(bp['whiskers'][1].get_data()[1][1],
           bp['whiskers'][3].get_data()[1][1])
minv = ax.get_ylim()[0]
ax.annotate('', xy=(0, maxv), xycoords = 'data',
            xytext=(1, maxv), textcoords = 'data',
            arrowprops = dict(arrowstyle = '-', ec = 'k',
                    connectionstyle = 'bar,fraction=0.2'))
_, pvalue = mannwhitneyu(h1_active, h1_inactive)
text_y1 = maxv+abs(maxv-minv)*0.12
ax.text(0.5, text_y1, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)
# pair 2
maxv = max(bp['whiskers'][5].get_data()[1][1],
           bp['whiskers'][7].get_data()[1][1])
minv = ax.get_ylim()[0]
ax.annotate('', xy=(2, maxv), xycoords = 'data',
            xytext=(3, maxv), textcoords = 'data',
            arrowprops = dict(arrowstyle = '-', ec = 'k',
                    connectionstyle = 'bar,fraction=0.2'))
_, pvalue = mannwhitneyu(hff_active, hff_inactive)
text_y2 = maxv+abs(maxv-minv)*0.19
ax.text(2.5, text_y2, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Number of linked enhancers', fontsize=6)
ax.set_xlabel('Genes in Lamina', fontsize=6)
ax.set_xticklabels(['Active', 'Inactive', 'Active', 'Inactive'],
                   rotation=45, ha='right', fontsize=5)
ax.spines['bottom'].set_linewidth(0.8)
ax.spines['left'].set_linewidth(0.8)
ax.xaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=0.8, labelsize=5, pad=2)
#ax.tick_params(axis='x', bottom=False, top=False, labelbottom=False, labeltop=False)
ax.set_axisbelow(True)
plt.savefig('enhancer_num-boxplot.svg', dpi=500, bbox_inches='tight')
plt.close()
