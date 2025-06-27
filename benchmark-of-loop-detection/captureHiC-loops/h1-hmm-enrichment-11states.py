import bisect, random, joblib, os, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LogNorm
from palettable.colorbrewer.diverging import RdYlBu_7_r

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def load_loops(infil):

    pchic_only = []
    shared = []
    with open(infil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2, label1, label2 = line.rstrip().split()
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            if label1 == 'anchor1':
                loop = (c2, s2, e2, c1, s1, e1)
            else:
                loop = (c1, s1, e1, c2, s2, e2)
            
            if label2 == '.':
                pchic_only.append(loop)
            else:
                shared.append(loop)
    
    return pchic_only, shared

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
        interval = abs(s2 - s1)
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

random_n = 1000
# load input data
hmm_fil = os.path.join('H1ESC.ChromHMM_8marks_12states.bed')
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
po_pchic_only, po_shared = load_loops('captureHiC-stratified-loops.po.bedpe')
pp_pchic_only, pp_shared = load_loops('captureHiC-stratified-loops.pp.bedpe')
reorg_loops = {
    'po_pchic_only': po_pchic_only,
    'po_shared': po_shared,
    'pp_pchic_only': pp_pchic_only,
    'pp_shared': pp_shared
}

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

state_briefs = ['AP', 'WP', 'PP', 'SE', 'WE', 'Insulator', 'TT', 'TE', 'WT', 'Repressed']
for i in sorted(reorg_loops):
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    farr = score_pool[i]
    farr[farr==0] = farr[farr>0].min()
    data = pd.DataFrame(farr, columns=state_briefs, index=state_briefs)
    cg = sns.heatmap(data, cmap='YlGnBu_r', annot=True, ax=ax, annot_kws={'size': 11},
         vmax=12, vmin=1, yticklabels=True, xticklabels=True, fmt='.1f')
    cbar = cg.figure.axes[-1]
    cbar.set_ylabel('Fold Change Enrichment', fontsize=16, rotation=270, labelpad=15)
    plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=16)
    plt.setp(cg.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=16)
    cg.set_clip_on(False)
    plt.savefig('{0}.hmm-fc.svg'.format(i), dpi=300, bbox_inches='tight')
    plt.close()