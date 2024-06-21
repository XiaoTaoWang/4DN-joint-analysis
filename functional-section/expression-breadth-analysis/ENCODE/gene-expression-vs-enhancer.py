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

def parse_exp(fil):

    F = open(fil, 'r')
    labels = F.readline().rstrip().split('\t')[5:]
    D = {}
    for line in F:
        parse = line.rstrip().split('\t')
        tmp = [float(i) for i in parse[5:]]
        gene_id = parse[4]
        D[gene_id] = {}
        for cell, tpm in zip(labels, tmp):
            D[gene_id][cell] = tpm
    
    F.close()

    return D

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

def read_rna_expression(exp_fil, gene_db, cell):

    exp = {}
    with open(exp_fil, 'r') as source:
        header = source.readline().rstrip().split('\t')
        cell_idx = header.index(cell)
        for line in source:
            parse = line.rstrip().split('\t')
            chrom, start, end, strand, ID, name = parse[:6]
            tpm = float(parse[cell_idx])
            if ID in gene_db:
                exp[ID] = tpm
    
    return exp

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
encode_fil = '/Users/xiaotaowang/workspace/ENCODE-ChIA-PET-2023/Gene-expressions/Gene-quants.matrix.log2.qn.avg.tsv'
gene_db = read_total_genes(exp_fil)
encode_exp_db = parse_exp(encode_fil)
queue = glob.glob('*.distal-enhancers.tsv')
# gene expressions
noenhancer = []
per10 = []
per20 = []
per30 = []
per40 = []
per50 = []
per60 = []
per70 = []
per80 = []
per90 = []
per100 = []
# number of tissues
noenhancer_num = []
per10_num = []
per20_num = []
per30_num = []
per40_num = []
per50_num = []
per60_num = []
per70_num = []
per80_num = []
per90_num = []
per100_num = []
for fil in queue:
    try:
        counts = read_enhancer_num(fil, gene_db)
        num_pool = np.r_[list(counts.values())]
        cell = os.path.split(fil)[1].split('.')[0]
        gene_exp = read_rna_expression(exp_fil, gene_db, cell)
    except:
        #print(cell, 'no RNA-Seq data available')
        continue
    
    if cell == 'H1':
        continue

    print(cell)
    # genes without distal enhancers
    for g in gene_exp:
        if not g in counts:
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            noenhancer.append(gene_exp[g])
            noenhancer_num.append(num_expressed)
    
    per10_cutoff = percentile(num_pool, 10, 0, 0)
    for g in counts:
        if (counts[g] > 0) and (counts[g] <= per10_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per10.append(gene_exp[g])
            per10_num.append(num_expressed)

    per20_cutoff = percentile(num_pool, 20, per10_cutoff, 10)
    for g in counts:
        if (counts[g] > per10_cutoff) and (counts[g] <= per20_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per20.append(gene_exp[g])
            per20_num.append(num_expressed)

    per30_cutoff = percentile(num_pool, 30, per20_cutoff, 20)
    for g in counts:
        if (counts[g] > per20_cutoff) and (counts[g] <= per30_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per30.append(gene_exp[g])
            per30_num.append(num_expressed)

    per40_cutoff = percentile(num_pool, 40, per30_cutoff, 30)
    for g in counts:
        if (counts[g] > per30_cutoff) and (counts[g] <= per40_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per40.append(gene_exp[g])
            per40_num.append(num_expressed)

    per50_cutoff = percentile(num_pool, 50, per40_cutoff, 40)
    for g in counts:
        if (counts[g] > per40_cutoff) and (counts[g] <= per50_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per50.append(gene_exp[g])
            per50_num.append(num_expressed)

    per60_cutoff = percentile(num_pool, 60, per50_cutoff, 50)
    for g in counts:
        if (counts[g] > per50_cutoff) and (counts[g] <= per60_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per60.append(gene_exp[g])
            per60_num.append(num_expressed)

    per70_cutoff = percentile(num_pool, 70, per60_cutoff, 60)
    for g in counts:
        if (counts[g] > per60_cutoff) and (counts[g] <= per70_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per70.append(gene_exp[g])
            per70_num.append(num_expressed)
    
    per80_cutoff = percentile(num_pool, 80, per50_cutoff, 70)
    for g in counts:
        if (counts[g] > per70_cutoff) and (counts[g] <= per80_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per80.append(gene_exp[g])
            per80_num.append(num_expressed)

    per90_cutoff = percentile(num_pool, 90, per50_cutoff, 80)
    for g in counts:
        if (counts[g] > per80_cutoff) and (counts[g] <= per90_cutoff):
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per90.append(gene_exp[g])
            per90_num.append(num_expressed)
    
    for g in counts:
        if counts[g] > per90_cutoff:
            num_expressed = len([encode_exp_db[g][_] for _ in encode_exp_db[g] if encode_exp_db[g][_]>3])
            per100.append(gene_exp[g])
            per100_num.append(num_expressed)

fig = plt.figure(figsize=(1.8, 1.5))
ax = fig.add_subplot(111)
bp = ax.boxplot([noenhancer, per10, per20, per30, per40, per50, per60, per70, per80, per90, per100],
                positions=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], sym='', widths=0.6,
                patch_artist=True, notch=False,
                whiskerprops={'linestyle':'--', 'linewidth':0.7})
print(list(map(len, [noenhancer, per10, per20, per30, per40, per50, per60, per70, per80, per90, per100])))
for patch in bp['boxes']:
    patch.set_facecolor('#FDBF6F')
    patch.set_linewidth(0.7)
for median in bp['medians']:
    median.set(color='k', linewidth=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('RNA-Seq log2(TPM + 1)', fontsize=6)
ax.set_xticklabels(['no\nenhancer', '(0, 10]', '(10, 20]', '(20, 30]', '(30, 40]',
                    '(40, 50]', '(50, 60]', '(60, 70]', '(70, 80]', '(80, 90]', '(90, 100]'],
                    fontsize=5, rotation=90, ha='center', va='top')
ax.set_xlabel('Percentile of the number of enhancers\nlinked to promoter', fontsize=6)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_axisbelow(True)
plt.savefig('exp-vs-number_of_enhancers.svg', dpi=500, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(1.8, 1.5))
ax = fig.add_subplot(111)
nbins = 50
sns.distplot(noenhancer_num, bins=nbins, kde=True, rug=False, color='b', ax=ax, hist=False, label='Without distal enhancers')
sns.distplot(per10_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(1/10), ax=ax, hist=False, label='(0, 10]')
sns.distplot(per20_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(2/10), ax=ax, hist=False, label='(10, 20]')
sns.distplot(per30_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(3/10), ax=ax, hist=False, label='(20, 30]')
sns.distplot(per40_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(4/10), ax=ax, hist=False, label='(30, 40]')
sns.distplot(per50_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(5/10), ax=ax, hist=False, label='(40, 50]')
sns.distplot(per60_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(6/10), ax=ax, hist=False, label='(50, 60]')
sns.distplot(per70_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(7/10), ax=ax, hist=False, label='(60, 70]')
sns.distplot(per80_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(8/10), ax=ax, hist=False, label='(70, 80]')
sns.distplot(per90_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(9/10), ax=ax, hist=False, label='(80, 90]')
sns.distplot(per100_num, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(10/10), ax=ax, hist=False, label='(90, 100]')
ax.set_xlim(-5, 125)
ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
ax.legend(frameon=False, fontsize=5)
ax.set_xlabel('Expression breadth (# of tissues)', fontsize=6)
ax.set_ylabel('Density', fontsize=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_axisbelow(True)
plt.savefig('expression-breadth.svg', dpi=500, bbox_inches='tight')
plt.close()
