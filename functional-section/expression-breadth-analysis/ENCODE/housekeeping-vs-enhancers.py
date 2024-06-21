import joblib, glob, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pyensembl import EnsemblRelease
import seaborn as sns
from matplotlib import ticker as mticker
from collections import Counter
import seaborn as sns

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def read_total_genes(exp_fil, min_tpm=1):
    
    genes = set()
    chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
    with open(exp_fil, 'r') as source:
        header = source.readline().rstrip().split('\t')
        for line in source:
            parse = line.rstrip().split('\t')
            chrom, start, end, strand, ID, name = parse[:6]
            if not chrom in chroms:
                continue
            tpm = max([float(t) for t in parse[6:]])
            if tpm <= min_tpm:
                continue
            
            genes.add(ID)
    
    return genes

def read_enhancer_num(fil, gene_db):

    genes = []
    with open(fil, 'r') as source:
        for line in source:
            gid, c, s, e, direction, loop_class = line.rstrip().split()
            if not gid in gene_db:
                continue
            genes.append(gid)
    
    counts = Counter(genes)

    return counts

def load_house_keeping():

    housekeeping = set()
    with open('human-housekeeping.txt', 'r') as source:
        for line in source:
            housekeeping.add(line.rstrip().split()[0])
    
    return housekeeping

def percentile(num_pool, per, low, low_per):

    v1 = np.percentile(num_pool, per)
    v2 = v1 - 1
    mask_1 = (num_pool > low) & (num_pool <= v1)
    mask_2 = (num_pool > low) & (num_pool <= v2)
    per_1 = mask_1.sum() / mask_1.size * 100
    per_2 = mask_2.sum() / mask_2.size * 100
    diff = per - low_per

    if abs(per_1 - diff) < abs(per_2 - diff):
        return v1
    else:
        return v2

exp_fil = '/Users/xiaotaowang/workspace/ENCODE-ChIA-PET-2023/Gene-expressions/Gene-quants.total-RNASeq.qn.avg.tsv'
gene_db = read_total_genes(exp_fil)
queue = glob.glob('*.distal-enhancers.tsv')
# collect genes with different number of distal enhancers
noenhancer = []
group1 = []
group2 = []
group3 = []
group4 = []
group5 = []
group6 = []
group7 = []
group8 = []
group9 = []
group10 = []
for fil in queue:
    counts = read_enhancer_num(fil, gene_db)
    num_pool = np.r_[list(counts.values())]
    for g in gene_db:
        if not g in counts:
            noenhancer.append(g)
    
    per10_cutoff = percentile(num_pool, 10, 0, 0)
    for g in counts:
        if counts[g] <= per10_cutoff:
            group1.append(g)

    for g in counts:
        if counts[g] > per10_cutoff:
            group2.append(g)
    
    per20_cutoff = percentile(num_pool, 20, per10_cutoff, 10)
    for g in counts:
        if counts[g] > per20_cutoff:
            group3.append(g)

    per30_cutoff = percentile(num_pool, 30, per20_cutoff, 20)
    for g in counts:
        if counts[g] > per30_cutoff:
            group4.append(g)

    per40_cutoff = percentile(num_pool, 40, per30_cutoff, 30)
    for g in counts:
        if counts[g] > per40_cutoff:
            group5.append(g)
    
    per50_cutoff = percentile(num_pool, 50, per40_cutoff, 40)
    for g in counts:
        if counts[g] > per50_cutoff:
            group6.append(g)
    
    per60_cutoff = percentile(num_pool, 60, per50_cutoff, 50)
    for g in counts:
        if counts[g] > per60_cutoff:
            group7.append(g)
    
    per70_cutoff = percentile(num_pool, 70, per60_cutoff, 60)
    for g in counts:
        if counts[g] > per70_cutoff:
            group8.append(g)
    
    per80_cutoff = percentile(num_pool, 80, per50_cutoff, 70)
    for g in counts:
        if counts[g] > per80_cutoff:
            group9.append(g)
    
    per90_cutoff = percentile(num_pool, 90, per50_cutoff, 80)
    for g in counts:
        if counts[g] > per90_cutoff:
            group10.append(g)

housekeeping = load_house_keeping()
background_avg = len(housekeeping) / len(gene_db)
group_pool = [noenhancer, group2, group3, group4,
              group5, group6, group7, group8, group9, group10]
sample_num_cutoffs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
M = np.zeros((len(sample_num_cutoffs), len(group_pool)))
for i in range(len(group_pool)):
    current_group = group_pool[i]
    counts = Counter(current_group)
    for j in sample_num_cutoffs:
        tmp_genes = set()
        for g in counts:
            if counts[g] > j:
                tmp_genes.add(g)
        
        tmp = housekeeping & tmp_genes
        ratio = len(tmp) / len(tmp_genes)
        fc = ratio / background_avg
        M[j, i] = fc

fig = plt.figure(figsize=(2, 2))
ax = fig.add_axes([0.1, 0.1, 0.7, 0.7])
sc = ax.pcolormesh(M, cmap='bwr', vmax=2.5, vmin=0.5)
ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
ax.set_xticklabels(['no\nenhancer', '>10', '>20', '>30', '>40',
                    '>50', '>60', '>70', '>80', '>90'], fontsize=5,
                    rotation=90, ha='center', va='top')
ax.set_xlabel('Enhancer number percentile', fontsize=6)
ax.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
ax.set_yticklabels(['>0', '>1', '>2', '>3', '>4', '>5', '>6', '>7', '>8', '>9', '>10'],
                   fontsize=5)
ax.set_ylabel('Number of biosamples', fontsize=6)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_axisbelow(True)
'''
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j+0.5, i+0.5, '{0:.2g}'.format(M[i,j]),
                ha='center', va='center', fontsize=4, color="w")
'''
ax = fig.add_axes([0.84, 0.65, 0.03, 0.15])
fig.colorbar(sc, cax=ax, ticks=[0.5, 2.5], format='%.3g')
ax.tick_params(labelsize=5)
ax.set_ylabel('Enrichment of\nhouse-keeping genes', fontsize=5, rotation=270, labelpad=15)

plt.savefig('housekeeping-vs-enhancers.svg', dpi=500, bbox_inches='tight')
plt.close()