import joblib, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def parse_ENCODE_exp(fil):

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

cell = 'HFFc6'
numenhancers_by_gid = {}
exp_by_gid = {}
numcells_by_gid = {}
exp_db = parse_ENCODE_exp('Gene-quants.matrix.log2.qn.avg.tsv')
fea_matrix, y, row_info, feature_labels = joblib.load('{0}-prediction-features.all.2k.pkl'.format(cell))
for i in range(len(row_info)):
    gid = row_info[i][4]
    if gid in exp_db:
        rna = np.log2(y[i][0]+1)
        num_enhancers = fea_matrix[i][5]
        num_expressed = len([exp_db[gid][_] for _ in exp_db[gid] if exp_db[gid][_]>3])
        numenhancers_by_gid[gid] = num_enhancers
        exp_by_gid[gid] = rna
        numcells_by_gid[gid] = num_expressed
        
print(len(exp_by_gid))

groups = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (6, 8),
    (8, 11),
    (11, np.inf)
]

exp_pool = []
numcells_pool = []
for l, u in groups:
    exp_tmp = []
    numcells_tmp = []
    for gid in exp_by_gid:
        if l <= numenhancers_by_gid[gid] < u:
            exp_tmp.append(exp_by_gid[gid])
            numcells_tmp.append(numcells_by_gid[gid])
    print(l, u, len(exp_tmp))
    exp_pool.append(exp_tmp)
    numcells_pool.append(numcells_tmp)

fig = plt.figure(figsize=(1.8, 1.2))
ax = fig.add_subplot(111)
bp = ax.boxplot(exp_pool, positions=[0, 1, 2, 3, 4, 5, 6, 7, 8], sym='', widths=0.6, patch_artist=True, notch=False,
                whiskerprops={'linestyle':'--', 'linewidth':0.7})
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
ax.set_xticklabels(['0', '1', '2', '3', '4', '5', '6-7', '8-10', '10+'], fontsize=5, ha='center')
ax.set_xlabel('Number of enhancers linked to promoter', fontsize=6)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=5, pad=2)
ax.set_axisbelow(True)
plt.savefig('{0}.number_of_enhancers.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()


colors = ['b'] + [plt.cm.Reds(i/8) for i in range(1, 9)]
labels = [
    '0',
    '1',
    '2',
    '3',
    '4',
    '5',
    '6-7',
    '8-10',
    '10+'
]
fig = plt.figure(figsize=(1.5, 1.2))
ax = fig.add_subplot(111)
nbins = 50
for numcells, color, label in zip(numcells_pool, colors, labels):
    sns.distplot(numcells, bins=nbins, kde=True, rug=False, color=color, ax=ax, hist=False, label=label)
ax.legend(frameon=False, fontsize=4)
ax.set_xlim(-5, 125)
ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
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
plt.savefig('{0}.expression-breadth.svg'.format(cell), dpi=500, bbox_inches='tight')
plt.close()