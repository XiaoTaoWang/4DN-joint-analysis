import bisect, random, glob, os, sys, joblib
import numpy as np

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

def get_loop_loci(loops):

    loci = set()
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        loci.add((c1, s1, e1))
        loci.add((c2, s2, e2))

    return loci

def generate_control_list(loops, chromfil, gapfil):

    gaps = parseGapFile(gapfil)
    chromsizes = parseChromLens(chromfil)

    random_list = []
    visited = set()
    for l in loops:
        c1, s1, e1, c2, s2, e2 = l
        res = max(e1-s1, e2-s2)
        check = True
        if c1 == c2:
            interval = s2 - s1
            while check:
                rp = random.randint(0, chromsizes[c1]-interval-res*2)
                tmp = [rp, rp+res]
                cache = check_in(tmp, gaps[c1])
                key = (c1, rp, rp+res, c1, rp+interval, rp+interval+res)
                if (len(cache)==0) and (not key in visited):
                    random_list.append(list(key))
                    visited.add(key)
                    check = False
        else:
            while check:
                rp1 = random.randint(0, chromsizes[c1]-res*2)
                rp2 = random.randint(0, chromsizes[c2]-res*2)
                tmp1 = [rp1, rp1+res]
                tmp2 = [rp2, rp2+res]
                cache1 = check_in(tmp1, gaps[c1])
                cache2 = check_in(tmp2, gaps[c2])
                key = (c1, rp1, rp1+res, c2, rp2, rp2+res)
                if (len(cache1)==0) and (len(cache2)==0) and (not key in visited):
                    random_list.append(list(key))
                    visited.add(key)
                    check = False
    
    return random_list

def overlap_with_chip(loop_loci, chip_peaks, min_len=15000):

    hit = 0
    for c, s, e in loop_loci:
        res = e - s
        if res < min_len:
            half = (min_len - res) // 2
            tmp = [max(0, s-half), e+half]
        else:
            tmp = [s, e]
        if not c in chip_peaks:
            continue
        cache = check_in(tmp, chip_peaks[c])
        if len(cache):
            hit += 1
    
    return hit

TF_fils = glob.glob(os.path.join('TF-peaks', '*.GRCh38.bed'))
chromfil = 'hg38.chrom.sizes'
gapfil = 'hg38.gap.txt'
random_n = 500
labels = joblib.load('consensus-clusters.new_labels.reassigned.pkl')
loops = joblib.load('H1ESC.umap-input.pkl')[0]
loops_by_clusters = {}
for i in sorted(set(labels)):
    if i == -1:
        continue
    idx = np.where(labels==i)[0]
    loops_by_clusters[i] = [list(loops[j][:6]) for j in idx]

score_pool = {}
for ID in sorted(loops_by_clusters):
    loops = loops_by_clusters[ID]
    loop_loci = list(get_loop_loci(loops))
    random_pool = []
    for i in range(random_n):
        random_loops = generate_control_list(loops, chromfil, gapfil)
        random_loci = list(get_loop_loci(random_loops))[:len(loop_loci)]
        random_pool.append(random_loci)
    
    Map = {}
    for f in TF_fils:
        path, fname = os.path.split(f)
        TF = fname.split('.')[1]
        chip_peaks = parseChIP(f)
        t_per = overlap_with_chip(loop_loci, chip_peaks) / len(loop_loci)
        r_per = []
        for i in range(random_n):
            total = len(random_pool[i])
            per = overlap_with_chip(random_pool[i], chip_peaks) / total
            r_per.append(per)
        enrich = t_per / np.mean(r_per)
        Map[TF] = [t_per, enrich]
        print('{0}, {1}, {2}, {3}'.format(ID, TF, t_per, enrich))
    
    score_pool[ID] = Map

joblib.dump(score_pool, 'H1-union-loops.TF-enrich.pkl')






