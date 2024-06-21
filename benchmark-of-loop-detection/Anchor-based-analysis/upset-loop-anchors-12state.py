import numpy as np
from upsetplot import plot
from upsetplot import from_memberships
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.colors import LogNorm
import seaborn as sns
import pandas as pd
import sys, bisect, random, os, joblib
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        if self.vmin == self.vmax:
            return np.ma.masked_array(np.interp(value, [self.vmin], [0.5]))
        
        if self.vmin < self.midpoint < self.vmax:
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        elif self.vmin >= self.midpoint:
            x, y = [self.vmin, self.vmax], [0.5, 1]
        elif self.vmax <= self.midpoint:
            x, y = [self.vmin, self.vmax], [0, 0.5]
            
        return np.ma.masked_array(np.interp(value, x, y))

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
            cache = overlap_peaks(tmp, gaps[c])
            key = (c, rp, rp+interval)
            if (len(cache)==0) and (not key in visited):
                random_list.append(list(key))
                visited.add(key)
                check = False
    
    return random_list

def overlap_peaks(p, List):

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

def check_in(p, List, mode='binary', min_len=15000):

    interval = p[1] - p[0]
    if interval < min_len:
        half = (min_len - interval) // 2
        tmp = [max(0, p[0]-half), p[1]+half]
        p = tmp
    interval = p[1] - p[0]

    cache = set()
    idx = max(0, bisect.bisect(List, p)-1)
    for q in List[idx:]:
        if q[1] <= p[0]:
            continue
        if q[0] >= p[1]:
            break
        minp, maxp = p[0], p[1]
        if q[0] > minp:
            minp = q[0]
        if q[1] < maxp:
            maxp = q[1]
        cache.add((minp, maxp))
    
    score = 0
    if mode=='binary':
        if len(cache):
            score = 1
    else:
        score = 0
        for s, e in cache:
            score += (e-s)
        score = score / interval

    
    return score

cell = sys.argv[1].split('.')[0]
anchors_by_platforms = {}
with open(sys.argv[1], 'r') as source:
    for line in source:
        parse = line.rstrip().split()
        tmp = parse[-1].split(',')
        labels = []
        for t in tmp:
            if t in ['MicroC-Feng', 'MicroC-Job']:
                labels.append('Micro-C')
            elif t in ['HiC-Feng', 'HiC-Job']:
                labels.append('Hi-C')
            else:
                labels.append(t)
        key = tuple(sorted(set(labels)))

        if not key in anchors_by_platforms:
            anchors_by_platforms[key] = []
        anchors_by_platforms[key].append((parse[0], int(parse[1]), int(parse[2])))

# data for upset plot
categories = [list(i) for i in anchors_by_platforms]
counts_by_class = [len(anchors_by_platforms[i]) for i in anchors_by_platforms]
print(counts_by_class)
data = from_memberships(
    categories,
    data = counts_by_class
)
fig = plt.figure(figsize=(4, 2))
plot(data, fig=fig, sort_categories_by='cardinality', sort_by=None)
plt.savefig('upset-anchors.svg', dpi=300, bbox_inches='tight')
plt.close()

# extract anchors in order
if cell == 'H1ESC':
    labels_in_order = ['Micro-C', 'Hi-C', 'ChIAPET-RNAPII', 'ChIAPET-CTCF', 'PLACSeq-H3K4me3'] # H1
else:
    labels_in_order = ['Micro-C', 'ChIAPET-RNAPII', 'Hi-C', 'ChIAPET-CTCF', 'PLACSeq-H3K4me3'] # HFF

key_in_order = []
for n in range(1, 6):
    queue = list(combinations(labels_in_order, n))
    for q in queue:
        key_in_order.append(q)

# data for fold enrichment
chromfil = 'hg38.chrom.sizes'
gapfil = 'hg38.gap.txt'
states_code = [
    '1_Active_Promoter',
    '2_Weak_Promoter',
    '3_Poised_Promoter',
    '4_Strong_Enhancer',
    '5_Strong_Enhancer',
    '6_Weak_Enhancer',
    '7_Insulator',
    '8_Transcriptional_Transition',
    '9_Transcriptional_Elongation',
    '10_Weak_Transcribed',
    '11_Polycomb_Repressed'
]
chromhmm_folder = '/Users/xtwang/workspace/4DN-Joint-Analysis/4DN-method-comparison/Loop-calls-HiC2.5/chromhmm-states-myrun/4DN_H1_HFF_8marks_12'
hmm_fil = os.path.join(chromhmm_folder, '{0}.ChromHMM_8marks_12states.bed'.format(cell))

states = parse_hmm(hmm_fil)
enrich_pool = []
for key in key_in_order:
    real_anchors = anchors_by_platforms[tuple(sorted(key))]
    score_pool = []
    for n in range(101):
        print(key, n)
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
                    scores.append(0)
                    continue
                scores.append(check_in([s, e], state[c], mode='binary', min_len=10000))
            arr[j] = np.mean(scores)
        score_pool.append(arr)
    
    score_pool = np.r_[score_pool]
    enrich_score = score_pool[0] / score_pool[1:].mean(axis=0)
    enrich_pool.append(enrich_score)

enrich_pool = np.r_[enrich_pool]
enrich_pool = enrich_pool.T

joblib.dump(enrich_pool, '{0}-anchors.hmm-enrich.12state.pkl'.format(cell))

fig = plt.figure(figsize=(15.6, 4.8))
ax = fig.add_subplot(111)
data = pd.DataFrame(enrich_pool, index=states_code)
cg = sns.heatmap(data, cmap='RdBu_r', annot=True, ax=ax, norm=MidpointNormalize(vmin=0, vmax=enrich_pool.max(), midpoint=1), annot_kws={'fontsize':7})
plt.setp(cg.yaxis.get_majorticklabels(), rotation=0, fontsize=10)
plt.setp(cg.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=10)
cg.set_clip_on(False)
plt.savefig('{0}.hmm-fc.12state.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()

                

