import sys, os, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def load_SAF_values(infil, per=25):

    by_gene = {}
    all_saf_values = []
    with open(infil, 'r') as source:
        source.readline()
        for line in source:
            ID, c, s, e, saf, gname, label = line.rstrip().split(',')
            all_saf_values.append(float(saf))
            by_gene[gname] = float(saf)
    
    low_cutoff = np.percentile(all_saf_values, per)
    high_cutoff = np.percentile(all_saf_values, 100-per)
    
    return by_gene, low_cutoff, high_cutoff

def load_enhancer_num(infil):

    loop_nums = {}
    with open(infil, 'r') as source:
        for line in source:
            ID, gname, num = line.rstrip().split()
            loop_nums[gname] = int(num)
    
    return loop_nums

SAF_folder = 'Frank'
cell = 'H1ESC'
colors = ['r', '#1F78B4', 'k']
loop_nums = load_enhancer_num('{0}.human-housekeeping.txt'.format(cell))
saf_fil = os.path.join(SAF_folder, '{0}_saf_gene.csv'.format(cell))
by_gene, low_cutoff, high_cutoff = load_SAF_values(saf_fil)
high_saf = []
low_saf = []
medium_saf = []
for g in loop_nums:
    if g in by_gene:
        if by_gene[g] > high_cutoff:
            high_saf.append(loop_nums[g])
        elif by_gene[g] < low_cutoff:
            low_saf.append(loop_nums[g])
        else:
            medium_saf.append(loop_nums[g])

print(len(high_saf), len(medium_saf), len(low_saf))

fig = plt.figure(figsize=(0.9, 1.3))
ax = fig.add_subplot(111)
bp = sns.violinplot(data=[high_saf, low_saf], ax=ax, palette=colors)
for violin in ax.collections[:]:
    violin.set_alpha(0.5)

print(mannwhitneyu(high_saf, low_saf))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Number of interacting enhancers', fontsize=6)
ax.set_xticklabels(['Class I', 'Class II'], rotation=45, ha='right', fontsize=5)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_axisbelow(True)
ax.set_xlabel('House-keeping genes', fontsize=6)
plt.savefig('{0}.enhancer-num-by-SAF.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()

