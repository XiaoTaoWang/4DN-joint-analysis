import bisect

def parse_hg38_ELS(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            if not c in D:
                D[c] = []
            D[c].append((s, e))
    
    for c in D:
        D[c].sort()
    
    return D

def parse_ChIP(fil):

    peaks = []
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            peaks.append((c, int(s), int(e)))
    
    return peaks

def bisect_search(p, ref):

    cache = set()
    if not p[0] in ref:
        return cache
    
    List = ref[p[0]]
    idx = max(0, bisect.bisect(List, p[1:])-1)
    for q in List[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add((p[0],)+q)
    
    return cache

ref_enhancers = parse_hg38_ELS('GRCh38-ELS.bed')
'''
cell = 'H1ESC'
chip_fils = [
    'H1ESC-ATACpeaks.4DNFI247OOFU.extended.bed'
]
'''
cell = 'HFFc6'
chip_fils = [
    'HFFc6-ATACpeaks.4DNFIWQJFZHS.extended.bed'
]

chip_peaks = []
for fil in chip_fils:
    chip_peaks.extend(parse_ChIP(fil))

pool = set()
for c, s, e in chip_peaks:
    cache = bisect_search((c, s, e), ref_enhancers)
    if len(cache):
        pool.add((c, s, e))

with open('{0}.enhancer-elements.bed'.format(cell), 'w') as out:
    for c, s, e in sorted(pool):
        out.write('\t'.join([c, str(s), str(e)])+'\n')
