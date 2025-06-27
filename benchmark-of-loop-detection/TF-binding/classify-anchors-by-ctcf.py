import bisect, joblib
import numpy as np

def get_loop_loci(loops):

    loci = set()
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        loci.add((c1, s1, e1))
        loci.add((c2, s2, e2))

    return loci

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

def parseChIP(chipfil):

    peaks = {}
    with open(chipfil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if '_' in chrom:
                continue
            if not chrom in peaks:
                peaks[chrom] = []
            sb = int(parse[1])
            eb = int(parse[2])
            peaks[chrom].append([sb, eb])
    
    for c in peaks:
        peaks[c].sort()
    
    return peaks

peaks = parseChIP('TF-peaks/H1ESC.CTCF.GRCh38.bed')
labels = joblib.load('consensus-clusters.new_labels.reassigned.pkl')
loops = joblib.load('H1ESC.umap-input.pkl')[0]
loop_loci_by_cluster = {}
for i in sorted(set(labels)):
    if i == -1:
        continue
    idx = np.where(labels==i)[0]
    loop_loci = get_loop_loci([list(loops[j][:6]) for j in idx])
    loop_loci_by_cluster[i] = loop_loci

with_ctcf = {}
without_ctcf = {}
min_len = 15000
for i in loop_loci_by_cluster:
    with_ctcf[i] = []
    without_ctcf[i] = []
    for l in loop_loci_by_cluster[i]:
        c, s, e = l
        res = e - s
        if res < min_len:
            half = (min_len - res) // 2
            tmp = [max(0, s-half), e+half]
        else:
            tmp = [s, e]

        if c in peaks:
            if len(check_in(tmp, peaks[c])):
                with_ctcf[i].append(l)
            else:
                without_ctcf[i].append(l)

joblib.dump(with_ctcf, 'with-CTCF.by_cluster.pkl')
joblib.dump(without_ctcf, 'without-CTCF.by_cluster.pkl')
