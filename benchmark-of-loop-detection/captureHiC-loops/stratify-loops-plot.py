import matplotlib
import numpy as np 
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.gridspec as gridspec

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def load_data(infil):

    collect = []
    with open(infil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            if parse[-1] == '.':
                label = 'PCHi-C'
            else:
                label = parse[-1]
            
            collect.append(label)
    
    counts = Counter(collect)

    return counts

colors = ['#7d423e', '#8d4945', '#9d4f4d', '#ae5859', '#bd6264', '#cb6e73', '#d48589', '#dc989c',
          '#896e3c', '#9d8243', '#b19749', '#c4a753', '#d5b25d', '#e1bc6b', '#e2c584', '#e4cc9a',
          '#6b7942', '#7b894d', '#8b9a59', '#9dad65', '#aebf70', '#bdce7e', '#c7d491', '#cfd9a2',
          '#3b3f77', '#454987', '#505498', '#5b5fa3', '#6569a9', '#7276b2', '#CCCCCC'
          ]
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
    'Micro-C', 'Hi-C,Micro-C', 'Hi-C', 'CTCF-ChIAPET,Hi-C,Micro-C',
    'CTCF-ChIAPET', 'CTCF-ChIAPET,Micro-C', 'CTCF-ChIAPET,Hi-C', 'PCHi-C'
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
    'CTCF & H3K4me3 & Hi-C & Pol2', 'CTCF & H3K4me3 & Hi-C',
    'Micro-C', 'Hi-C & Micro-C', 'Hi-C', 'CTCF & Hi-C & Micro-C',
    'CTCF', 'CTCF & Micro-C', 'CTCF & Hi-C', 'PCHi-C only'
]

counts_po = load_data('captureHiC-stratified-loops.po.bedpe')
counts_pp = load_data('captureHiC-stratified-loops.pp.bedpe')
arr_po = []
arr_pp = []
for k in keys:
    arr_po.append(counts_po[k]/sum(counts_po.values()))
    arr_pp.append(counts_pp[k]/sum(counts_pp.values()))

data = np.array([arr_po, arr_pp])
data_cum = data.cumsum(axis=1)

fig = plt.figure(figsize=(3, 2.2))
gs = gridspec.GridSpec(1, 6, wspace=0.35, top=0.9, bottom=0.1, left=0.1,
                       right=0.9)
ax = fig.add_subplot(gs[0:2])
handles = []
legend_labels = []
labels = ['PO', 'PP']
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
ax.set_ylabel('Percentage of PCHi-C interactions', fontsize=7)
plt.savefig('stratify.svg', dpi=300, bbox_inches='tight')
plt.close()