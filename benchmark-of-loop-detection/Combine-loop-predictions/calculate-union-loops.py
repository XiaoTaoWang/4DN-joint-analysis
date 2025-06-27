from collections import defaultdict
import numpy as np

def parse_bedpe(fil):

    pos1 = defaultdict(list)
    pos2 = defaultdict(list)
    loops = defaultdict(list)
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if chrom in ['chrM']:
                continue

            c1, s1, e1, c2, s2, e2 = parse[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            p1 = (s1 + e1) // 2
            p2 = (s2 + e2) // 2
            pos1[chrom].append(p1)
            pos2[chrom].append(p2)
            loops[chrom].append((c1, s1, e1, c2, s2, e2))
    
    for c in pos1:
        pos1[c] = np.r_[pos1[c]]
        pos2[c] = np.r_[pos2[c]]
    
    return pos1, pos2, loops

def check_in(p1, p2, ref1, ref2, ref_loops, c, thre=10000):

    hit = []
    if not c in ref1:
        return hit
    
    tmp1 = ref1[c]
    tmp2 = ref2[c]
    thre = min(thre, abs(p2-p1)*0.2)
    dis1 = np.abs(p1 - tmp1)
    dis2 = np.abs(p2 - tmp2)
    mask = (dis1 <= thre) & (dis2 <= thre)
    idx = np.where(mask)[0]
    if idx.size > 0:
        hit = [ref_loops[c][i] for i in idx]
    
    return hit

### H1ESC
fil_list = [
    'H1ESC_CTCF.Feng.loops.peak_based.bedpe',
    'H1ESC_RNAPII.Feng.loops.peak_based.bedpe',
    'H1ESC.PLAC-H3K4me3.peakachu.clustered.peak_based.bedpe',
    'H1ESC-MicroC.Feng.loops.2kb.bedpe',
    'H1ESC-HiC.Feng.loops.2kb.bedpe',
    'H1ESC_CTCF.Feng.loops.2kb.bedpe',
    'H1ESC_RNAPII.Feng.loops.2kb.bedpe',
    'H1ESC-MicroC.Job.loops.4DNFI3RMWQ85.5kb.bedpe',
    'H1ESC-HiC.Job.loops.4DNFI4PMUZNU.5kb.bedpe',
    'H1ESC-MicroC.Feng.loops.5kb.bedpe',
    'H1ESC-HiC.Feng.loops.5kb.bedpe',
    'H1ESC_CTCF.Feng.loops.5kb.bedpe',
    'H1ESC_RNAPII.Feng.loops.5kb.bedpe',
    'H1ESC.PLAC-H3K4me3.peakachu.clustered.5kb.bedpe',
    'H1ESC-MicroC.Job.loops.4DNFI3RMWQ85.10kb.bedpe',
    'H1ESC-HiC.Job.loops.4DNFI4PMUZNU.10kb.bedpe',
    'H1ESC-MicroC.Feng.loops.10kb.bedpe',
    'H1ESC-HiC.Feng.loops.10kb.bedpe',
    'H1ESC_PLACSeq.Bing.loops.10kb.bedpe'
]
names = ['CTCF-ChIAPET', 'RNAPII-ChIAPET', 'H3K4me3-PLACSeq', 'Micro-C', 'Hi-C',
         'CTCF-ChIAPET', 'RNAPII-ChIAPET', 'Micro-C', 'Hi-C', 'Micro-C', 'Hi-C',
         'CTCF-ChIAPET', 'RNAPII-ChIAPET', 'H3K4me3-PLACSeq', 'Micro-C', 'Hi-C', 'Micro-C',
         'Hi-C', 'H3K4me3-PLACSeq']
outfil = 'H1ESC.union-loops.bedpe'
'''
### HFFc6
fil_list = [
    'HFFc6_CTCF.Feng.loops.peak_based.bedpe',
    'HFFc6_RNAPII.Feng.loops.peak_based.bedpe',
    'HFFc6.PLAC-H3K4me3.peakachu.clustered.peak_based.bedpe',
    'HFFc6-MicroC.Feng.loops.2kb.bedpe',
    'HFFc6-HiC.Feng.loops.2kb.bedpe',
    'HFFc6_CTCF.Feng.loops.2kb.bedpe',
    'HFFc6_RNAPII.Feng.loops.2kb.bedpe',
    'HFFc6-MicroC.Job.loops.4DNFIIQP46FO.5kb.bedpe',
    'HFFc6-HiC.Job.loops.4DNFI212GCZW.5kb.bedpe',
    'HFFc6-MicroC.Feng.loops.5kb.bedpe',
    'HFFc6-HiC.Feng.loops.5kb.bedpe',
    'HFFc6_CTCF.Feng.loops.5kb.bedpe',
    'HFFc6_RNAPII.Feng.loops.5kb.bedpe',
    'HFFc6.PLAC-H3K4me3.peakachu.clustered.5kb.bedpe',
    'HFFc6-MicroC.Job.loops.4DNFIIQP46FO.10kb.bedpe',
    'HFFc6-HiC.Job.loops.4DNFI212GCZW.10kb.bedpe',
    'HFFc6-MicroC.Feng.loops.10kb.bedpe',
    'HFFc6-HiC.Feng.loops.10kb.bedpe',
    'HFFc6_PLACSeq.Bing.loops.10kb.bedpe'
]
names = ['CTCF-ChIAPET', 'RNAPII-ChIAPET', 'H3K4me3-PLACSeq', 'Micro-C', 'Hi-C',
         'CTCF-ChIAPET', 'RNAPII-ChIAPET', 'Micro-C', 'Hi-C', 'Micro-C', 'Hi-C',
         'CTCF-ChIAPET', 'RNAPII-ChIAPET', 'H3K4me3-PLACSeq', 'Micro-C', 'Hi-C', 'Micro-C',
         'Hi-C', 'H3K4me3-PLACSeq']
outfil = 'HFFc6.union-loops.bedpe'
'''
loop_sets = [parse_bedpe(f) for f in fil_list]
thre = 15000

cache = set()
union_loops = []
for i in range(len(loop_sets)):
    query1, query2, query_loops = loop_sets[i]
    for c in query1:
        for p1, p2, tmp_loop in zip(query1[c], query2[c], query_loops[c]):
            if tmp_loop in cache: # current loop have been added in previous steps
                continue

            ID = [0] * len(fil_list)
            ID[i] = 1
            for j in range(len(loop_sets)):
                if i==j:
                    continue
                ref1, ref2, ref_loops = loop_sets[j]
                hit = check_in(p1, p2, ref1, ref2, ref_loops, c, thre=thre)
                if len(hit):
                    ID[j] = 1
                for h_loop in hit:
                    cache.add(h_loop)

            cache.add(tmp_loop)
            
            loop_label = ','.join(sorted(set([names[ii] for ii in range(len(ID)) if ID[ii]==1])))
            if sum(ID) >= 1:
                union_loops.append(tmp_loop + (loop_label,))

with open(outfil, 'w') as out:
    for loop in sorted(union_loops):
        out.write('\t'.join(list(map(str, loop)))+'\n')

    
    
