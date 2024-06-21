import joblib, cooler, os
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from multiprocess import Pool

def _expected_core(pool_args):
    
    hic_pool, c, maxdis, balance = pool_args

    expected = {} # average over each genomic distance
    hic = hic_pool.matrix(balance=balance, sparse=True).fetch(c)
    
    #tmp = np.array(hic.sum(axis=0)).ravel() > 0 # filter out gap regions
    weights = hic_pool.bins().fetch(c)['weight'].values
    tmp = np.isfinite(weights) & (weights > 0) # filter out gap regions
    n = hic.shape[0]
    maxdis = min(n-1, maxdis)
    # Assign values for each genomic distance
    for i in range(1, maxdis+1):
        valid = tmp[:-i] * tmp[i:]
        current = hic.diagonal(i)[valid]
        current = current[np.isfinite(current)]
        if current.size > 0:
            expected[i] = [current.sum(), current.size]
    
    return expected

def calculate_expected(uri, chroms, maxdis=5000000, balance=True, nproc=1):

    hic_pool = cooler.Cooler(uri)
    res = hic_pool.binsize
    maxdis = maxdis // res
    args = []
    for c in chroms:
        args.append((hic_pool, c, maxdis, balance))
    
    # Allocate processes
    if nproc==1:
        results = list(map(_expected_core, args))
    else:
        pool = Pool(nproc)
        results = pool.map(_expected_core, args)
        pool.close()
        pool.join()

    expected = {}
    for i in range(1, maxdis+1):
        nume = 0
        denom = 0
        for extract in results:
            if i in extract:
                nume += extract[i][0]
                denom += extract[i][1]
        if nume > 0:
            expected[i] = nume / denom
    
    ini_x = sorted(expected)
    ini_exp = [expected[i] for i in ini_x]
    IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
    IR.fit(ini_x, ini_exp)
    x = list(range(1, maxdis+1))
    expected = dict(zip(x, IR.predict(x)))
    
    return expected

'''
queue = [
    ('Hi-C', 'HIC-H1ESC-DpnII-allReps.4DNFI82R42AD.mcool'),
    ('Micro-C', 'MicroC-H1ESC-MNase-allReps.mcool'),
    ('H3K4me3 PLAC-Seq', 'PLAC-H1ESC-H3K4me3-allReps.mcool'),
    ('CTCF ChIA-PET', 'ChIAPET-H1ESC-CTCF-allReps.mcool'),
    ('RNAPII ChIA-PET', 'ChIAPET-H1ESC-RNAPII-allReps.mcool'),
    ('SPRITE', 'H1ESC.2_1000.two_over_n.mcool'),
    ('SPRITE_2_10', 'H1ESC.2_10.two_over_n.mcool')
]
'''
queue = [
    ('GAM', 'H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at25Kb.mcool')
]

h1_loops, hff_loops, labels = joblib.load('umap-6clusters-cache.pkl')
labels = labels[:len(h1_loops)]
loops = h1_loops
label_map = {-1:-1, 0:-1, 1:-1, 2:0, 3:-1, 4:1, 5:2, 6:3, 7:4, 8:5}
new_labels = []
for i in labels:
    new_labels.append(label_map[i])
labels = np.r_[new_labels]

scores_pool = np.zeros((len(queue), 6))
res = 25000
min_size = 30000
max_size = 3000000
n_process = 6
chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
platforms = []
for i in range(len(queue)):
    platform, mcool = queue[i]
    platforms.append(platform)
    print(platform)
    if platform in ['Hi-C', 'Micro-C', 'SPRITE', 'SPRITE_2_10']:
        balance = True
    else:
        balance = False
    uri = '{0}::resolutions/{1}'.format(mcool, res)
    clr = cooler.Cooler(uri)
    expected = calculate_expected(uri, chroms, maxdis=max_size, balance=balance, nproc=n_process)
    for j in range(6):
        print('cluster: {0}'.format(j+1))
        class_label = j
        idx = np.where(labels==class_label)[0]
        current_loops = {}
        for k in idx:
            c1, s1, e1, c2, s2, e2 = loops[k][:6]
            if (s2 - s1 > max_size) or (s2 - s1 < min_size):
                continue
            if not c1 in current_loops:
                current_loops[c1] = []
            current_loops[c1].append((s1, e1, s2, e2))
        
        scores = []
        for c in current_loops:
            M = clr.matrix(balance=balance, sparse=True).fetch(c).tocsr()
            for s1, e1, s2, e2 in current_loops[c]:
                bins1 = range(s1//res, (e1+res-1)//res)
                bins2 = range(s2//res, (e2+res-1)//res)
                maxv = 0
                binpair = None
                for b1 in bins1:
                    for b2 in bins2:
                        if M[b1, b2] > maxv:
                            maxv = M[b1, b2]
                            binpair = (b1, b2)

                if not binpair is None:
                    dis = abs(binpair[1] - binpair[0])
                    if dis in expected:
                        normed_value = maxv / expected[dis]
                        scores.append(normed_value)

        scores_pool[i, j] = np.mean(scores)

cluster_labels = [str(i) for i in range(1, 7)]
data = pd.DataFrame(scores_pool, columns=cluster_labels, index=platforms)
joblib.dump(data, 'H1.contact-strength.GAM.pkl')