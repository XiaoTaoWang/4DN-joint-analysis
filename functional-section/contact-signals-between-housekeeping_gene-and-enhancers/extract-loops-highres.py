import cooler, sys, pickle
from multiprocess import Pool
import numpy as np
from sklearn.isotonic import IsotonicRegression

def parse_real_loops(fil):

    loops = {}
    with open(fil, 'r') as source:
        for i, line in enumerate(source):
            c1, s1, e1, c2, s2, e2, label = line.rstrip().split()
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            if not c1 in loops:
                loops[c1] = []
            loops[c1].append((s1, e1, s2, e2, label, i))
    
    return loops

def _expected_core(pool_args):
    
    hic_pool, c, maxdis, balance = pool_args

    expected = {} # average over each genomic distance
    hic = hic_pool.matrix(balance=balance, sparse=True).fetch(c).tocsr()
    
    #tmp = np.array(hic.sum(axis=0)).ravel() > 0 # filter out gap regions
    weights = hic_pool.bins().fetch(c)['weight'].values
    tmp = np.isfinite(weights) & (weights > 0) # filter out gap regions
    n = hic.shape[0]
    maxdis = min(n-1, maxdis)
    # Assign values for each genomic distance
    for i in range(maxdis+1):
        if i == 0:
            valid = tmp
        else:
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
    for i in range(maxdis+1):
        nume = 0
        denom = 0
        for extract in results:
            if i in extract:
                nume += extract[i][0]
                denom += extract[i][1]
        if nume > 0:
            expected[i] = nume / denom
    
    IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
    ini_x = np.r_[sorted(expected)]
    ini_exp = np.r_[[expected[i] for i in ini_x]]
    _d = np.where(ini_exp > 0)[0]
    IR.fit(ini_x[_d], ini_exp[_d])
    x = list(range(maxdis+1))
    expected = dict(zip(x, IR.predict(x)))
    
    return expected

def distance_normaize_core(sub, x_arr, y_arr, expected):

    D = y_arr - x_arr[:,np.newaxis]
    D = np.abs(D)
    min_dis = D.min()
    max_dis = D.max()
    if max_dis > max(expected):
        return
    else:
        exp_sub = np.zeros(sub.shape)
        for d in range(min_dis, max_dis+1):
            xi, yi = np.where(D==d)
            exp_sub[xi, yi] = expected[d]
        normed = sub / exp_sub

        return normed

def distance_normalize(uri, loops, expected, width, max_apart=2000000, min_size=200000,
    balance=False):

    lib = cooler.Cooler(uri)
    chroms = loops.keys()

    res = lib.binsize
    square_n = (2*width + 1) ** 2
    arrs = {}
    for c in sorted(chroms):
        if not c in lib.chromnames:
            continue

        M = lib.matrix(balance=balance, sparse=True).fetch(c).tocsr()
        for s1, e1, s2, e2, label, idx in loops[c]:
            if s2 - s1 < min_size:
                continue
            if s2 - s1 > max_apart:
                continue
            
            bins1 = range(s1//res, (e1+res-1)//res)
            bins2 = range(s2//res, (e2+res-1)//res)
            maxv = -1
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
                        if label == 'downstream':
                            arrs[idx] = arr[width] # centered at enhancer loci
                        else:
                            arrs[idx] = arr[:,width]
    
    return arrs

uri = sys.argv[1]
cell = sys.argv[2]
outpickle = sys.argv[3]
'''
balance = sys.argv[3]

if balance == 'ICE':
    balance = True
else:
    balance = False
'''
balance = False
chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
width = 15
clr = cooler.Cooler(uri)
max_apart = 1000000
min_size = 10*clr.binsize
n_process = 4

real_fil = '{0}.EP-loops.housekeeping.bedpe'.format(cell)
real_loops = parse_real_loops(real_fil)
expected = calculate_expected(uri, chroms, maxdis=max_apart, balance=balance, nproc=n_process)
print(expected[0], expected[1], expected[2])
arrs_real = distance_normalize(uri, real_loops, expected, width, max_apart=max_apart, min_size=min_size,
                               balance=balance)

with open(outpickle, 'wb') as out:
    pickle.dump(arrs_real, out)
