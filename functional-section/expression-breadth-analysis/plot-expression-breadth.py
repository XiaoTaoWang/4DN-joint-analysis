import matplotlib
import numpy as np
import matplotlib.pyplot as plt
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

cell = 'HFFc6'
infil = '{0}.distal-enhancers.tsv'.format(cell)
outfil = '{0}.exp-breadth.svg'.format(cell)
encode_exp_db = parse_exp('Gene-quants.matrix.log2.qn.avg.tsv')
num_exp_1 = []
num_exp_2 = []
num_exp_3 = []
num_exp_4 = []
num_exp_5 = []
num_exp_6 = []
num_exp_7 = []
num_exp_8 = []
num_exp_9 = []
with open(infil, 'r') as source:
    source.readline()
    for line in source:
        parse = line.strip().split('\t')
        if parse[0] in ['chrM', 'chrY']:
            continue
        if parse[6] == 'no':
            continue

        gid = parse[3]
        tpm = float(parse[4])
        enhancers = set()
        for e in parse[7].split(';'):
            enhancers.add(tuple(e.split(',')[1:]))
        n = len(enhancers)
        num_expressed = len([encode_exp_db[gid][_] for _ in encode_exp_db[gid] if encode_exp_db[gid][_]>3])
        if parse[7] == '.':
            num_exp_1.append(num_expressed)
        elif n == 1:
            num_exp_2.append(num_expressed)
        elif n == 2:
            num_exp_3.append(num_expressed)
        elif n == 3:
            num_exp_4.append(num_expressed)
        elif n == 4:
            num_exp_5.append(num_expressed)
        elif n == 5:
            num_exp_6.append(num_expressed)
        elif n in [6, 7]:
            num_exp_7.append(num_expressed)
        elif n in [8, 9, 10]:
            num_exp_8.append(num_expressed)
        elif n > 8:
            num_exp_9.append(num_expressed)

fig = plt.figure(figsize=(1.6, 1.2))
ax = fig.add_subplot(111)
nbins = 50
sns.distplot(num_exp_1, bins=nbins, kde=True, rug=False, color='b', ax=ax, hist=False, label='Without distal enhancers')
sns.distplot(num_exp_2, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(1/8), ax=ax, hist=False, label='1 distal enhancer')
sns.distplot(num_exp_3, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(2/8), ax=ax, hist=False, label='2 distal enhancers')
sns.distplot(num_exp_4, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(3/8), ax=ax, hist=False, label='3 distal enhancers')
sns.distplot(num_exp_5, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(4/8), ax=ax, hist=False, label='4 distal enhancer')
sns.distplot(num_exp_6, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(5/8), ax=ax, hist=False, label='5 distal enhancers')
sns.distplot(num_exp_7, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(6/8), ax=ax, hist=False, label='6-7 distal enhancers')
sns.distplot(num_exp_8, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(7/8), ax=ax, hist=False, label='8-10 distal enhancers')
sns.distplot(num_exp_9, bins=nbins, kde=True, rug=False, color=plt.cm.Reds(8/8), ax=ax, hist=False, label='10+ distal enhancers')
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
