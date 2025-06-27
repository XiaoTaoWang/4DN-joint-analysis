import bisect, random, os, joblib, matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.diverging import RdYlBu_7_r

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

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

def generate_control_list(true_loci, chromfil, gapfil):

    gaps = parseGapFile(gapfil)
    chromsizes = parseChromLens(chromfil)

    random_list = []
    visited = set()
    for c, s, e in true_loci:
        interval = e - s
        check = True
        while check:
            rp = random.randint(0, chromsizes[c]-interval-1)
            tmp = [rp, rp+interval]
            cache = check_in(tmp, gaps[c])
            key = (c, rp, rp+interval)
            if (len(cache)==0) and (not key in visited):
                random_list.append(list(key))
                visited.add(key)
                check = False
    
    return random_list

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

hmm_fil = os.path.join('H1ESC.ChromHMM_8marks_12states.bed')
states = parse_hmm(hmm_fil)
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
state_labels = [
    'Active Promoter',
    'Weak Promoter',
    'Poised Promoter',
    'Strong Enhancer',
    'Weak Enhancer',
    'Insulator',
    'Transcriptional Transition',
    'Transcriptional Elongation',
    'Weakly Transcribed',
    'Polycomb Repressed'
]
random_n = 501
loop_loci_by_clusters = joblib.load('with-CTCF.by_cluster.pkl')
outfil = 'with-CTCF.ChromHMM-enrich.pkl'
'''
enrich_pool = []
for ID in sorted(loop_loci_by_clusters):
    real_anchors = loop_loci_by_clusters[ID]
    score_pool = []
    for n in range(random_n):
        print(ID, n)
        if n == 0:
            anchors = real_anchors
        else:
            anchors = generate_control_list(real_anchors, chromfil, gapfil)
        
        arr = np.zeros(len(states_code))
        for j, s_n in enumerate(states_code):
            state = states[s_n]
            scores = []
            for c, s, e in anchors:
                if not c in state:
                    continue
                check = overlap_hmm([s, e], state[c])
                if check:
                    scores.append(1)
                else:
                    scores.append(0)
            arr[j] = np.mean(scores)
        
        score_pool.append(arr)
    
    score_pool = np.r_[score_pool]
    enrich_score = score_pool[0] / score_pool[1:].mean(axis=0)
    enrich_pool.append(enrich_score)

enrich_pool = np.r_[enrich_pool]
joblib.dump(enrich_pool, outfil)
'''
column_labels = ['AP', 'WP', 'PP', 'SE', 'WE', 'Insulator', 'TT', 'TE', 'WT', 'Repressed']
enrich_pool = joblib.load(outfil)
normed = []
for arr in enrich_pool:
    arr = arr / np.mean(arr)
    normed.append(arr)
normed = np.r_[normed]
fig = plt.figure(figsize=(3, 1))
ax = fig.add_subplot(111)
data = pd.DataFrame(normed[2:], index=['C3', 'C4', 'C5', 'C6'],
                    columns=column_labels)
cg = sns.heatmap(data, cmap=RdYlBu_7_r.get_mpl_colormap(), annot=enrich_pool[2:],
                 ax=ax, annot_kws={'fontsize':5}, vmax=3)
plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=6)
plt.setp(cg.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=6)
cg.set_clip_on(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_ylabel('Loop Clusters', fontsize=6)
ax.set_xlabel('ChromHMM States', fontsize=6)
ax.set_title('CTCF-bound loops', fontsize=6)
ax.xaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
ax.yaxis.set_tick_params(width=0.5, length=1.5, pad=1, labelsize=6)
plt.savefig(outfil.replace('.pkl', '.svg'), dpi=300, bbox_inches='tight')
plt.close()