import bisect, random, joblib, os, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LogNorm
from collections import defaultdict
from palettable.colorbrewer.diverging import RdYlBu_7_r
'''
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)
'''
def parseChromLens(chromfil):
    
    chromsizes = {}
    with open(chromfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if '_' in chrom:
                continue
            chromsizes[chrom] = int(parse[1])
            
    return chromsizes

def parseGapFile(gapfil):
    
    gaps = {}
    with open(gapfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[1]
            if '_' in chrom:
                continue
            if not chrom in gaps:
                gaps[chrom] = []
            sb = int(parse[2])
            eb = int(parse[3])
            gaps[chrom].append([sb, eb+1])

    for c in gaps:
        gaps[c].sort()
    
    return gaps

def check_in(p, List):

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        cache.add(tuple(q))
    
    return cache

def generate_control_list(loops, chromfil, gapfil):

    gaps = parseGapFile(gapfil)
    chromsizes = parseChromLens(chromfil)

    random_list = []
    visited = set()
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        res = e1 - s1
        interval = s2 - s1
        check = True
        while check:
            rp = random.randint(0, chromsizes[c1]-interval-res*2)
            tmp = [rp, rp+interval]
            cache = check_in(tmp, gaps[c1])
            key = (c1, rp, rp+res, c1, rp+interval, rp+interval+res)
            if (len(cache)==0) and (not key in visited):
                random_list.append(list(key))
                visited.add(key)
                check = False
    
    return random_list

def parse_hmm(fil):

    hmm = {}
    with open(fil, 'r') as source:
        for line in source:
            if line.startswith('#'):
                continue
            parse = line.rstrip().split('\t')
            chrom, s, e, label = parse[:4]
            label = '_'.join(label.split('_')[1:])
            s, e = int(s), int(e)
            if not label in hmm:
                hmm[label] = {}
            if not chrom in hmm[label]:
                hmm[label][chrom] = []
            hmm[label][chrom].append([s, e])
    
    for label in hmm:
        for chrom in hmm[label]:
            hmm[label][chrom].sort()
    
    return hmm

def overlap_hmm(p, List, min_len=10000):

    interval = p[1] - p[0]
    if interval < min_len:
        half = (min_len - interval) // 2
        tmp = [max(0, p[0]-half), p[1]+half]
        p = tmp

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        cache.add(tuple(q))
    
    if len(cache):
        return True
    else:
        return False

cell = 'HFFc6'
random_n = 500
# load input data
hmm_fil = os.path.join('{0}.ChromHMM_8marks_12states.bed'.format(cell))
chromfil = 'hg38.chrom.sizes'
gapfil = 'hg38.gap.txt'
states_code = [
    'Active_Promoter',
    'Weak_Promoter',
    'Poised_Promoter',
    'Strong_Enhancer',
    'Weak_Enhancer',
    'Insulator',
    'Transcriptional_Transition',
    'Transcriptional_Elongation',
    'Weak_Transcribed',
    'Polycomb_Repressed'
]
states = parse_hmm(hmm_fil)
loops = joblib.load('{0}.umap-input.pkl'.format(cell))[0]
labels = joblib.load('consensus-clusters.new_labels.pkl')
reorg_loops = {}
for i in sorted(set(labels)):
    if i == -1:
        continue
    idx = np.where(labels==i)[0]
    reorg_loops[i] = [loops[j][:6] for j in idx]
'''
score_pool = {}
for c_i in reorg_loops:
    real_loops = reorg_loops[c_i]
    tmp_pool = []
    for n in range(random_n+1):
        print(c_i, n)
        if n == 0:
            loops = real_loops
        else:
            loops = generate_control_list(real_loops, chromfil, gapfil)
        
        arr = np.zeros((len(states_code), len(states_code)))
        for c1, s1, e1, c2, s2, e2 in loops:
            for i, r_s in enumerate(states_code):
                state1 = states[r_s]
                if not c1 in state1:
                    continue
                check1 = overlap_hmm([s1, e1], state1[c1])
                if not check1:
                    continue
                for j, c_s in enumerate(states_code):
                    state2 = states[c_s]
                    if c2 in state2:
                        check2 = overlap_hmm([s2, e2], state2[c2])
                        if check2:
                            arr[i, j] += 1

        tmp_pool.append(arr)

    tmp_pool = np.r_[tmp_pool]
    true_scores = tmp_pool[0]
    random_scores = tmp_pool[1:]
    farr = np.ones_like(true_scores)
    for i in range(true_scores.shape[0]):
        for j in range(true_scores.shape[1]):
            t_s = true_scores[i, j]
            r_s = random_scores[:, i, j]
            u = np.mean(r_s)
            farr[i, j] = t_s / u
    
    score_pool[c_i] = farr

joblib.dump(score_pool, '{0}.hmm-intermediate.by_umap_class.pkl'.format(cell))
'''
score_pool = joblib.load('{0}.hmm-intermediate.by_umap_class.pkl'.format(cell))
state_briefs = ['AP', 'WP', 'PP', 'SE', 'WE', 'Insulator', 'TT', 'TE', 'WT', 'Repressed']
cluster_labels = list(map(str, sorted(reorg_loops)))
top_n = 3
collect_pair_labels = []
for i in sorted(reorg_loops):
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    farr = score_pool[i]
    farr[farr==0] = farr[farr>0].min()
    data = pd.DataFrame(farr, columns=states_code, index=states_code)
    #cg = sns.heatmap(data, cmap='YlGnBu_r', annot=True, ax=ax, norm=LogNorm(vmin=farr.min(), vmax=farr.max()))
    cg = sns.heatmap(data, cmap='YlGnBu_r', annot=True, ax=ax, norm=LogNorm(vmin=farr.min(), vmax=farr.max()),
         vmax=10, yticklabels=True, xticklabels=True, fmt='.1f')
    cbar = cg.figure.axes[-1]
    cbar.set_ylabel('Fold Change Enrichment', fontsize=16, rotation=270, labelpad=15)
    plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=16)
    plt.setp(cg.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=16)
    cg.set_clip_on(False)
    plt.savefig('{0}.{1}.hmm-fc.svg'.format(cell, i), dpi=300, bbox_inches='tight')
    plt.close()

    # plot top n features
    sort_table = []
    for ii in range(len(farr)):
        for jj in range(len(farr)):
            sort_table.append((farr[ii,jj], (ii,jj)))
    sort_table.sort(reverse=True)
    for ii in range(top_n):
        tmp_pairs = tuple([state_briefs[sort_table[ii][1][0]], state_briefs[sort_table[ii][1][1]]])
        tmp_pairs = tuple(sorted(tmp_pairs))
        tmp_label = '-'.join(tmp_pairs)
        if not tmp_label in collect_pair_labels:
            collect_pair_labels.append(tmp_label)

collect_pair_labels = [
    'AP-AP', 'AP-WP', 'AP-TE', 'TE-TT', 'AP-SE', 'SE-WP', 'SE-SE',
    'PP-PP', 'PP-Repressed', 'PP-WP', 'Insulator-Insulator', 'Insulator-WP',
    'Insulator-PP', 'AP-TT', 'TT-TT'
]
fc_matrix = np.zeros((len(collect_pair_labels), len(reorg_loops)))
for i in sorted(reorg_loops):
    farr = score_pool[i]
    farr[farr==0] = farr[farr>0].min()
    sort_table = []
    for ii in range(len(farr)):
        for jj in range(len(farr)):
            tmp_pairs = tuple([state_briefs[ii], state_briefs[jj]])
            tmp_pairs = tuple(sorted(tmp_pairs))
            tmp_label = '-'.join(tmp_pairs)
            if tmp_label in collect_pair_labels:
                idx = collect_pair_labels.index(tmp_label)
                fc_matrix[idx, i] = max(farr[ii,jj], farr[jj,ii])

fig = plt.figure(figsize=(1.8, 2.1))
ax = fig.add_subplot(111)
fc_norm = fc_matrix / fc_matrix.max(axis=0)
data = pd.DataFrame(fc_norm, columns=cluster_labels, index=collect_pair_labels)
cg = sns.heatmap(data, cmap=RdYlBu_7_r.get_mpl_colormap(), annot=fc_matrix, ax=ax,
                 annot_kws={'fontsize':5}, fmt='.1f',
                 yticklabels=True, xticklabels=True)
plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=6)
plt.setp(cg.xaxis.get_majorticklabels(), rotation=0, fontsize=6)
#cg.xaxis.set_ticks_position('top')
#cg.yaxis.set_ticks_position('right')
cg.set_clip_on(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(0.5)
#ax.spines['left'].set_linewidth(0.5)
ax.xaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
ax.yaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
plt.savefig('{0}.clusters-features.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()