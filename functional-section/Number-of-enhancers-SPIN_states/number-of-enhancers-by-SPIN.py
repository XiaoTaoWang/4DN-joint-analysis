import joblib, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def read_SPIN_states(rnaseq_fil, cell='H1ESC'):

    D = {}
    with open(rnaseq_fil, 'r') as source:
        source.readline()
        for line in source:
            parse = line.rstrip().split('\t')
            gid = parse[0].split('.')[0]
            h1_state = parse[6]
            hff_state = parse[9]
            if cell == 'H1ESC':
                D[gid] = h1_state
            else:
                D[gid] = hff_state
    
    return D

def read_num_enhancers(cell):

    fea_matrix, y, row_info, feature_labels = joblib.load('{0}-prediction-features.all.2k.pkl'.format(cell))
    D = {}
    for i in range(len(row_info)):
        gid = row_info[i][4]
        num = int(fea_matrix[i][5])
        D[gid] = num
    
    return D

cell = 'H1ESC'
num_by_gene = read_num_enhancers(cell)
SPIN_by_gene = read_SPIN_states('RNA-seq_TPM_FPKM.v29.added_loop_info.tsv', cell)
names = ['Speckle', 'Interior_Act1', 'Interior_Act2', 'Interior_Act3',
        'Interior_Repr1', 'Interior_Repr2', 'Near_Lm1', 'Near_Lm2', 'Lamina']
num_by_states = defaultdict(list)
for gid in num_by_gene:
    if gid in SPIN_by_gene:
        key = SPIN_by_gene[gid]
        num_by_states[key].append(num_by_gene[gid])

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
