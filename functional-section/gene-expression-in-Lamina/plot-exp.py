import joblib, bisect, matplotlib
import numpy as np
from collections import defaultdict, Counter
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon
from pyensembl import EnsemblRelease

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

hg38_db = EnsemblRelease(94, species='human')
hg38_db.download()
hg38_db.index()

h1_total, h1_enhancers = joblib.load('H1ESC.lamina-genes.pkl')
hff_total, hff_enhancers = joblib.load('HFFc6.lamina-genes.pkl')
h1_no_enhancers = h1_total - h1_enhancers
hff_no_enhancers = hff_total - hff_enhancers
h1_tpm_enhancers = np.r_[[g[5] for g in h1_enhancers]]
hff_tpm_enhancers = np.r_[[g[5] for g in hff_enhancers]]
h1_tpm_no = np.r_[[g[5] for g in h1_no_enhancers]]
hff_tpm_no = np.r_[[g[5] for g in hff_no_enhancers]]
h1_tpm_enhancers = np.log2(h1_tpm_enhancers+1)
hff_tpm_enhancers = np.log2(hff_tpm_enhancers+1)
h1_tpm_no = np.log2(h1_tpm_no+1)
hff_tpm_no = np.log2(hff_tpm_no+1)

colors = ['#7FC97F', '#999999', '#BEAED4', '#999999']
fig = plt.figure(figsize=(1.4, 1.2))
ax = fig.add_subplot(111)
bp = ax.boxplot([h1_tpm_enhancers, h1_tpm_no, hff_tpm_enhancers, hff_tpm_no],
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
_, pvalue = mannwhitneyu(h1_tpm_enhancers, h1_tpm_no)
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
_, pvalue = mannwhitneyu(hff_tpm_enhancers, hff_tpm_no)
text_y2 = maxv+abs(maxv-minv)*0.19
ax.text(2.5, text_y2, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('RNA-Seq log2(TPM+1)', fontsize=6)
ax.set_xlabel('Genes in Lamina', fontsize=6)
ax.set_xticklabels(['with enhancers', 'no enhancers', 'with enhancers', 'no enhancers'],
                   rotation=45, ha='right', fontsize=5)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
#ax.tick_params(axis='x', bottom=False, top=False, labelbottom=False, labeltop=False)
ax.set_axisbelow(True)
plt.savefig('tpm-boxplot.svg', dpi=500, bbox_inches='tight')
plt.close()

# output gene list
h1_gene_list = set()
with open('H1-gene-list.txt', 'w') as out:
    for g in h1_enhancers:
        out.write(hg38_db.gene_name_of_gene_id(g[4])+'\n')
        h1_gene_list.add(hg38_db.gene_name_of_gene_id(g[4]))

hff_gene_list = set()
with open('HFF-gene-list.txt', 'w') as out:
    for g in hff_enhancers:
        out.write(hg38_db.gene_name_of_gene_id(g[4])+'\n')
        hff_gene_list.add(hg38_db.gene_name_of_gene_id(g[4]))