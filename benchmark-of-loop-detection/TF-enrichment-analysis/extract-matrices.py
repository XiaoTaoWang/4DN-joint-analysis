import numpy as np
import joblib, pyBigWig, sys

def parse_loops(fil):

    loops = {}
    with open(fil, 'r') as source:
        source.readline()
        for line in source:
            parse = line.rstrip().split()
            c1, s1, e1, c2, s2, e2, labels, ID = parse
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            ID = int(ID)
            if not ID in loops:
                loops[ID] = []
            loops[ID].append((c1, s1, e1, c2, s2, e2))
    
    return loops

#loops_by_clusters = parse_loops('H1ESC.union-loops.with-cluster-info.bedpe')
#bw = pyBigWig.open('H1-CTCF.ENCFF332TNJ.bigWig')
loops_by_clusters = parse_loops(sys.argv[1])
bw = pyBigWig.open(sys.argv[2])

width = 50000
nbin_w = 50
nbin_p = 5
matrices_by_clusters = {}
for ID in loops_by_clusters:
    if ID == -1:
        continue
    matrices_1 = []
    matrices_2 = []
    for c1, s1, e1, c2, s2, e2 in loops_by_clusters[ID]:
        if s1 - width <= 0:
            continue
        if e1 + width >= bw.chroms()[c1]:
            continue
        if s2 - width <= 0:
            continue
        if e2 + width >= bw.chroms()[c2]:
            continue
        
        # anchor 1
        part1 = bw.stats(c1, s1-width, s1, nBins=nbin_w)
        part2 = bw.stats(c1, s1, e1, nBins=nbin_p)
        part3 = bw.stats(c1, e1, e1+width, nBins=nbin_w)
        anchor1 = []
        for t in part1+part2+part3:
            if t is None:
                anchor1.append(0)
            else:
                anchor1.append(t)
        anchor1 = np.r_[anchor1]
        matrices_1.append(anchor1)

        # anchor 2
        part1 = bw.stats(c2, s2-width, s2, nBins=nbin_w)
        part2 = bw.stats(c2, s2, e2, nBins=nbin_p)
        part3 = bw.stats(c2, e2, e2+width, nBins=nbin_w)
        anchor2 = []
        for t in part1+part2+part3:
            if t is None:
                anchor2.append(0)
            else:
                anchor2.append(t)
        anchor2 = np.r_[anchor2]
        matrices_2.append(anchor2)
    
    matrices_1 = np.r_[matrices_1]
    matrices_2 = np.r_[matrices_2]
    matrices_by_clusters[ID] = [matrices_1, matrices_2]
'''
for ID in matrices_by_clusters:
    select_idx = []
    matrices = matrices_by_clusters[ID]
    for i in range(len(matrices[0])):
        check_sum = 0
        for j in range(len(matrices)):
            tmp = matrices[j][i]
            if tmp.sum() > 0:
                check_sum += 1
        if check_sum == len(matrices):
            select_idx.append(i)
    
    new_matrices = []
    for m in matrices:
        new_matrices.append(m[select_idx])
    
    print('cluster {0}: total {1}, after filtering {2}'.format(ID, len(matrices[0]), len(new_matrices[0])))
    matrices_by_clusters[ID] = new_matrices
'''
joblib.dump(matrices_by_clusters, sys.argv[3], compress=('xz', 3))

            

