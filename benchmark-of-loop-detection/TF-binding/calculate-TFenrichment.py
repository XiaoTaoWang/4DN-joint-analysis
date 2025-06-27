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

def generate_control_list(loop_loci, chromfil, gapfil):

    gaps = parseGapFile(gapfil)
    chromsizes = parseChromLens(chromfil)
    random_list = []
    visited = set()
    for c, s, e in loop_loci:
        res = e - s
        check = True
        while check:
            rp = random.randint(0, chromsizes[c]-res*2)
            tmp = [rp, rp+res]
            cache = check_in(tmp, gaps[c])
            key = (c, rp, rp+res)
            if (len(cache)==0) and (not key in visited):
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
loop_loci_by_clusters = joblib.load('without-CTCF.by_cluster.pkl')
outfil = 'without-CTCF.TF-enrich.pkl'

score_pool = {}
for ID in sorted(loop_loci_by_clusters):
    loop_loci = loop_loci_by_clusters[ID]
    random_pool = []
    for i in range(random_n):
        random_loci = generate_control_list(loop_loci, chromfil, gapfil)
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

joblib.dump(score_pool, outfil)