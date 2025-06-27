import matplotlib, joblib
import numpy as np 
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.gridspec as gridspec

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def extract_specific_cluster(loops, cluster_labels, ID):

    collect = []
    for loop, ci in zip(loops, cluster_labels):
        label = ','.join(loop[-1])
        if ci == ID:
            collect.append(label)
    
    counts = Counter(collect)

    return counts

colors = ['#7d423e', '#8d4945', '#9d4f4d', '#ae5859', '#bd6264', '#cb6e73', '#d48589', '#dc989c',
          '#896e3c', '#9d8243', '#b19749', '#c4a753', '#d5b25d', '#e1bc6b', '#e2c584', '#e4cc9a',
          '#6b7942', '#7b894d', '#8b9a59', '#9dad65', '#aebf70', '#bdce7e', '#c7d491', '#cfd9a2',
          '#3b3f77', '#454987', '#505498', '#5b5fa3', '#6569a9', '#7276b2', '#898cbe']
keys = [
    'H3K4me3-PLACSeq', 'H3K4me3-PLACSeq,RNAPII-ChIAPET', 'RNAPII-ChIAPET',
    'H3K4me3-PLACSeq,Hi-C,Micro-C,RNAPII-ChIAPET', 'CTCF-ChIAPET,H3K4me3-PLACSeq,Hi-C,Micro-C,RNAPII-ChIAPET',
    'H3K4me3-PLACSeq,Hi-C,Micro-C', 'CTCF-ChIAPET,H3K4me3-PLACSeq,Hi-C,Micro-C',
    'H3K4me3-PLACSeq,Micro-C,RNAPII-ChIAPET', 'H3K4me3-PLACSeq,Micro-C', 'H3K4me3-PLACSeq,Hi-C',
    'H3K4me3-PLACSeq,Hi-C,RNAPII-ChIAPET', 'Micro-C,RNAPII-ChIAPET', 'CTCF-ChIAPET,Hi-C,Micro-C,RNAPII-ChIAPET',
    'CTCF-ChIAPET,H3K4me3-PLACSeq,Micro-C,RNAPII-ChIAPET', 'CTCF-ChIAPET,H3K4me3-PLACSeq,Micro-C',
    'CTCF-ChIAPET,H3K4me3-PLACSeq', 'CTCF-ChIAPET,Micro-C,RNAPII-ChIAPET', 'CTCF-ChIAPET,H3K4me3-PLACSeq,RNAPII-ChIAPET',
    'CTCF-ChIAPET,RNAPII-ChIAPET', 'Hi-C,RNAPII-ChIAPET', 'Hi-C,Micro-C,RNAPII-ChIAPET',
    'CTCF-ChIAPET,H3K4me3-PLACSeq,Hi-C,RNAPII-ChIAPET', 'CTCF-ChIAPET,H3K4me3-PLACSeq,Hi-C',
    'CTCF-ChIAPET,Hi-C,RNAPII-ChIAPET',
    'Micro-C', 'Hi-C,Micro-C', 'Hi-C', 'CTCF-ChIAPET,Hi-C,Micro-C',
    'CTCF-ChIAPET', 'CTCF-ChIAPET,Micro-C', 'CTCF-ChIAPET,Hi-C'
]
category_names = [
    'H3K4me3', 'H3K4me3 & Pol2', 'Pol2',
    'H3K4me3 & Hi-C & Micro-C & Pol2', 'CTCF & H3K4me3 & Hi-C & Micro-C & Pol2',
    'H3K4me3 & Hi-C & Micro-C', 'CTCF & H3K4me3 & Hi-C & Micro-C',
    'H3K4me3 & Micro-C & Pol2', 'H3K4me3 & Micro-C', 'H3K4me3 & Hi-C',
    'H3K4me3 & Hi-C & Pol2', 'Micro-C & Pol2', 'CTCF & Hi-C & Micro-C & Pol2',
    'CTCF & H3K4me3 & Micro-C & Pol2', 'CTCF & H3K4me3 & Micro-C',
    'CTCF & H3K4me3', 'CTCF & Micro-C & Pol2', 'CTCF & H3K4me3 & Pol2',
    'CTCF & Pol2', 'Hi-C & Pol2', 'Hi-C & Micro-C & Pol2',
    'CTCF & H3K4me3 & Hi-C & Pol2', 'CTCF & H3K4me3 & Hi-C', 'CTCF & Hi-C & Pol2',
    'Micro-C', 'Hi-C & Micro-C', 'Hi-C', 'CTCF & Hi-C & Micro-C',
    'CTCF', 'CTCF & Micro-C', 'CTCF & Hi-C'
]

cell = 'HFFc6'
cluster_labels = joblib.load('consensus-clusters.new_labels.reassigned.pkl')
loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]
counts_1 = extract_specific_cluster(loops, cluster_labels, ID=0)
counts_2 = extract_specific_cluster(loops, cluster_labels, ID=1)
counts_3 = extract_specific_cluster(loops, cluster_labels, ID=2)
counts_4 = extract_specific_cluster(loops, cluster_labels, ID=3)
counts_5 = extract_specific_cluster(loops, cluster_labels, ID=4)
counts_6 = extract_specific_cluster(loops, cluster_labels, ID=5)
ratios_1 = []
ratios_2 = []
ratios_3 = []
ratios_4 = []
ratios_5 = []
ratios_6 = []
for k in keys:
    ratios_1.append(counts_1[k]/sum(counts_1.values()))
    ratios_2.append(counts_2[k]/sum(counts_2.values()))
    ratios_3.append(counts_3[k]/sum(counts_3.values()))
    ratios_4.append(counts_4[k]/sum(counts_4.values()))
    ratios_5.append(counts_5[k]/sum(counts_5.values()))
    ratios_6.append(counts_6[k]/sum(counts_6.values()))

data = np.array([ratios_1, ratios_2, ratios_3, ratios_4, ratios_5, ratios_6])
data_cum = data.cumsum(axis=1)

fig = plt.figure(figsize=(3, 2.2))
gs = gridspec.GridSpec(1, 6, wspace=0.35, top=0.9, bottom=0.1, left=0.1,
                       right=0.9)
ax = fig.add_subplot(gs[0:2])
handles = []
legend_labels = []
labels = ['1', '2', '3', '4', '5', '6']
for i, (name, color) in enumerate(zip(category_names, colors)):
    heights = data[:, i]
    starts = data_cum[:, i] - heights
    b = ax.bar(labels, heights, bottom=starts, width=0.77, color=color)
    handles.append(b)
    legend_labels.append(name)

ax.set_xticklabels(labels)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)

ax.xaxis.set_tick_params(width=1, labelsize=7, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=7, pad=2)

ax.legend(handles[::-1], legend_labels[::-1], ncol=2, bbox_to_anchor=(1.05, -0.02, 1, 1), loc='lower left',
          labelspacing=0.3, frameon=True, fontsize=6)
ax.set_xlabel('Loop Cluster ID', fontsize=7)
ax.set_ylabel('Percentage', fontsize=7)
plt.savefig('stratify.svg', dpi=300, bbox_inches='tight')
plt.close()