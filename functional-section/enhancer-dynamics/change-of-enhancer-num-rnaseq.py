from pyensembl import EnsemblRelease
import bisect, joblib
from collections import defaultdict, Counter
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.stats import mannwhitneyu, wilcoxon
from scipy import stats

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def read_data(fil):

    num_by_ID = {}
    rna_by_ID = {}
    fea_matrix, y, row_info, feature_labels = joblib.load(fil)
    for i in range(len(row_info)):
        gid = row_info[i][4]
        if (y[i][0] > 0) and (fea_matrix[i][5] > 0):
            num_by_ID[gid] = fea_matrix[i][5]
            rna_by_ID[gid] = y[i][0]
    
    return num_by_ID, rna_by_ID


h1_num, h1_rna = read_data('H1ESC-prediction-features.all.2k.pkl')
hff_num, hff_rna = read_data('HFFc6-prediction-features.all.2k.pkl')

group1 = [] # H1>HFF
group2 = [] # H1=HFF
group3 = [] # H1<HFF
for ID in h1_rna:
    if ID in hff_rna:
        fc = h1_rna[ID] / hff_rna[ID]
        if h1_num[ID] > hff_num[ID]:
            group1.append(np.log2(fc))
        elif h1_num[ID] < hff_num[ID]:
            group3.append(np.log2(fc))
        else:
            group2.append(np.log2(fc))

group1 = np.r_[group1]
group2 = np.r_[group2]
group3 = np.r_[group3]
colors = ['#7FC97F', '#BEAED4', '#FDC086']
fig = plt.figure(figsize=(1.4, 1.3))
ax = fig.add_subplot(111)
bp = ax.boxplot([group1, group2, group3], positions=[0, 1, 2], sym='',
                widths=0.6, patch_artist=True, notch=False,
                whiskerprops={'linestyle':'--', 'linewidth':0.7})
maxv = max(bp['whiskers'][1].get_data()[1][1],
           bp['whiskers'][3].get_data()[1][1],
           bp['whiskers'][5].get_data()[1][1])
minv = ax.get_ylim()[0]
# pair 1
ax.annotate('', xy=(0, maxv), xycoords = 'data',
            xytext=(1, maxv), textcoords = 'data',
            arrowprops = dict(arrowstyle = '-', ec = 'k',
                    connectionstyle = 'bar,fraction=0.2'))
_, pvalue = mannwhitneyu(group1, group2)
text_y1 = maxv+abs(maxv-minv)*0.21
ax.text(0.5, text_y1, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)
# pair 2
maxv = maxv+abs(maxv-minv)*0.24
ax.annotate('', xy=(1, maxv), xycoords = 'data',
            xytext=(2, maxv), textcoords = 'data',
            arrowprops = dict(arrowstyle = '-', ec = 'k',
                    connectionstyle = 'bar,fraction=0.2'))
_, pvalue = mannwhitneyu(group2, group3)
text_y2 = maxv+abs(maxv-minv)*0.17
ax.text(1.5, text_y2, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)
'''
# pair 3
maxv = maxv+abs(maxv-minv)*0.2
ax.annotate('', xy=(0, maxv), xycoords = 'data',
            xytext=(2, maxv), textcoords = 'data',
            arrowprops = dict(arrowstyle = '-', ec = 'k',
                    connectionstyle = 'bar,fraction=0.1'))
_, pvalue = mannwhitneyu(group1, group3)
text_y3 = maxv+abs(maxv-minv)*0.14
ax.text(1, text_y3, 'P={0:.3g}'.format(pvalue), ha = 'center', va = 'center', fontsize=5)
ax.set_ylim(minv, text_y3)
'''

ax.set_ylim(minv, 8)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_linewidth(0.7)
for median in bp['medians']:
    median.set(color='k', linewidth=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('RNA-Seq log2(H1/HFF)', fontsize=6)
ax.set_xticklabels(['H1>HFF', 'H1=HFF', 'H1<HFF'], fontsize=6, ha='right', rotation=45)
ax.set_xlabel('Number of linked enhancers', fontsize=6)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.set_axisbelow(True)
plt.savefig('fc_by_num_groups.rna-seq.svg', dpi=500, bbox_inches='tight')
plt.close()