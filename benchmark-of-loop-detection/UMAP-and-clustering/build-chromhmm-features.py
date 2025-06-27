import numpy as np
import sys, bisect, os, joblib

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

def parse_bedpe(fil):

    loops = []
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if chrom in ['chrM', 'chrY']:
                continue
            c1, s1, e1, c2, s2, e2 = parse[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            key = tuple(parse[-1].split(','))
            loops.append((c1, s1, e1, c2, s2, e2, key))
    
    return loops

def check_in(p, List, mode='binary', min_len=10000):

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
cell = 'HFFc6'
hmm_fil = '{0}.ChromHMM_8marks_12states.bed'.format(cell)
states = parse_hmm(hmm_fil)
loops = parse_bedpe('{0}.union-loops.bedpe'.format(cell))
M = []
filtered_loops = []
for c1, s1, e1, c2, s2, e2, key in loops:
    scores = []
    for s_n in states_code:
        state = states[s_n]
        if not c1 in state:
            scores.append(0)
        else:
            scores.append(check_in([s1, e1], state[c1], mode='fraction', min_len=10000))
    
    for s_n in states_code:
        state = states[s_n]
        if not c2 in state:
            scores.append(0)
        else:
            scores.append(check_in([s2, e2], state[c2], mode='fraction', min_len=10000))
    
    if (max(scores[:len(states_code)]) > 0) and (max(scores[len(states_code):]) > 0):
        M.append(scores)
        filtered_loops.append((c1, s1, e1, c2, s2, e2, key))

M = np.r_[M]
normM = (M - M.mean(axis=0)) / M.std(axis=0)

joblib.dump([filtered_loops, normM], '{0}.umap-input.pkl'.format(cell), compress=('xz', 3))