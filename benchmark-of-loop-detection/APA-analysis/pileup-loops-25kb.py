import cooler, sys, pickle
from multiprocess import Pool
import numpy as np
from sklearn.isotonic import IsotonicRegression

def parse_loops(fil):

    loops = {}
    with open(fil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2 = line.rstrip().split()[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            if not c1 in loops:
                loops[c1] = []
            loops[c1].append((s1, e1, s2, e2))
    
    return loops

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

def distance_normaize_core(sub, x_arr, y_arr, expected):

    D = y_arr - x_arr[:,np.newaxis]
    D = np.abs(D)
    min_dis = D.min()
    max_dis = D.max()
    if (max_dis > max(expected)) or (min_dis == 0):
        return
    else:
        exp_sub = np.zeros(sub.shape)
        for d in range(min_dis, max_dis+1):
            xi, yi = np.where(D==d)
            exp_sub[xi, yi] = expected[d]
        normed = sub / exp_sub

        return normed


def distance_normalize(uri, loopfil, expected, max_apart=2000000, min_size=200000,
    balance=False):

    lib = cooler.Cooler(uri)
    loops = parse_loops(loopfil)
    chroms = loops.keys()

    res = lib.binsize
    width = 10
    square_n = (2*width + 1) ** 2
    arrs = []
    for c in chroms:
        if not c in lib.chromnames:
            continue

        M = lib.matrix(balance=balance, sparse=True).fetch(c).tocsr()
        for s1, e1, s2, e2 in loops[c]:
            if s2 - s1 < min_size:
                continue
            if s2 - s1 > max_apart:
                continue
            
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
                b1, b2 = binpair
                arr = M[(b1-width):(b1+width+1), (b2-width):(b2+width+1)].toarray()
                if arr.size == square_n:
                    arr[np.isnan(arr)] = 0
                    arr = arr.astype(np.float64)
                    _x = np.arange(b1-width, b1+width+1)
                    _y = np.arange(b2-width, b2+width+1)
                    arr = distance_normaize_core(arr, _x, _y, expected)
                    if not arr is None:
                        arrs.append(arr)
    
    return arrs

def pileup(arrs):

    mean_arr = np.r_[[np.mean(arr) for arr in arrs]]
    p99 = np.percentile(mean_arr, 99)
    p1 = np.percentile(mean_arr, 1)
    mask = (mean_arr < p99) & (mean_arr > p1)
    pool = np.r_[arrs]
    avg = pool[mask].mean(axis=0)

    return avg

uri = sys.argv[1]
loop_fil = sys.argv[2]
outpickle = sys.argv[3]
platform = sys.argv[4]

if platform in ['Hi-C', 'Micro-C', 'SPRITE']:
    balance = True
else:
    balance = False

chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
max_apart = 5000000
min_size = 250000
n_process = 5
expected = calculate_expected(uri, chroms, maxdis=max_apart, balance=balance, nproc=n_process)
arrs = distance_normalize(uri, loop_fil, expected, max_apart=max_apart, min_size=min_size, balance=balance)
#avg = pileup(arrs)

with open(outpickle, 'wb') as out:
    pickle.dump(arrs, out)
