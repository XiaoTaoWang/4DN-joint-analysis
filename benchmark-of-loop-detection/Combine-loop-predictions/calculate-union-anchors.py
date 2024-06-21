from collections import Counter
import numpy as np

def parse_bedpe(fil):

    anchors = {}
    with open(fil, 'r') as source:
        if 'Job' in fil:
            source.readline()
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0]
            if chrom in ['chrM']:
                continue
            c1, s1, e1, c2, s2, e2 = parse[:6]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            if not chrom in anchors:
                anchors[chrom] = set()
            anchors[chrom].add((c1, s1, e1))
            anchors[chrom].add((c2, s2, e2))
    
    for c in anchors:
        anchors[c] = sorted(anchors[c])

    pos = {}
    for c in anchors:
        pos[c] = np.r_[[(a[1]+a[2])//2 for a in anchors[c]]]
    
    return anchors, pos

def check_in(p, ref_pos, ref_anchors, c, thre=10000):

    hit = []
    if not c in ref_anchors:
        return hit
    tmp = ref_pos[c]
    dis = np.abs(p - tmp)
    mask = (dis <= thre)
    idx = np.where(mask)[0]
    if idx.size > 0:
        hit = [ref_anchors[c][i] for i in idx]
    
    return hit
'''
### H1ESC
fil_list = [
    'H1ESC_RNAPII.Feng.loops.bedpe',
    'H1ESC_CTCF.Feng.loops.bedpe',
    'H1ESC-MicroC.Feng.loops.bedpe',
    'H1ESC-HiC.Feng.loops.bedpe',
    'H1ESC-MicroC.Job.loops.4DNFI3RMWQ85.bedpe',
    'H1ESC-HiC.Job.loops.4DNFI4PMUZNU.bedpe',
    'H1ESC-PLAC-H3K4me3.conserved.bedpe'
]
outfil = 'H1ESC.union-anchors.bed'
'''
### HFFc6
fil_list = [
    'HFFc6_RNAPII.Feng.loops.bedpe',
    'HFFc6-CTCF.Feng.loops.bedpe',
    'HFFc6-MicroC.Feng.loops.bedpe',
    'HFFc6-HiC.Feng.loops.bedpe',
    'HFFc6-MicroC.Job.loops.4DNFIIQP46FO.bedpe',
    'HFFc6-HiC.Job.loops.4DNFI212GCZW.bedpe',
    'HFFc6-PLAC-H3K4me3.conserved.bedpe'
]
outfil = 'HFFc6.union-anchors.bed'

names = ['ChIAPET-RNAPII', 'ChIAPET-CTCF', 'MicroC-Feng', 'HiC-Feng', 'MicroC-Job', 'HiC-Job', 'PLACSeq-H3K4me3']

anchors = [parse_bedpe(f) for f in fil_list]
thre = 15000

cache = set()
union_anchors = []
for i in range(len(anchors)):
    query_anchors, query_pos = anchors[i]
    for c in query_anchors:
        for coord, p in zip(query_anchors[c], query_pos[c]):
            if coord in cache: # current loop have been added in previous steps
                continue
            ID = [0, 0, 0, 0, 0, 0, 0]
            ID[i] = 1
            for j in range(len(anchors)):
                if i==j:
                    continue
                ref_anchors, ref_pos = anchors[j]
                hit = check_in(p, ref_pos, ref_anchors, c, thre=thre)
                if len(hit):
                    ID[j] = 1
                for tmp in hit:
                    cache.add(tmp)

            cache.add(coord)
            
            loop_label = ','.join([names[ii] for ii in range(len(ID)) if ID[ii]==1])
            if sum(ID) >= 1:
                union_anchors.append(coord + (loop_label,))

with open(outfil, 'w') as out:
    for anchor in sorted(union_anchors):
        out.write('\t'.join(list(map(str, anchor)))+'\n')

    
    
