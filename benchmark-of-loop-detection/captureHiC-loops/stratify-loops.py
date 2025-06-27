import matplotlib
import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, euclidean

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def parse_loops(fil):

    D = {}
    info_map = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            c1, s1, e1, c2, s2, e2, label = parse[:7]
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            p1 = (s1 + e1) // 2
            p2 = (s2 + e2) // 2
            if not c1 in D:
                D[c1] = []
            D[c1].append((s1, e1, s2, e2, p1, p2, e2-s1))
            info_map[(c1, s1, e1, c2, s2, e2)] = label
    
    for c1 in D:
        D[c1] = np.r_[D[c1]]

    return D, info_map

loop_fil1 = 'GSE86189.significant-pp.H1ESC.hg38.bedpe'
loop_fil2 = 'H1ESC.union-loops.bedpe'
outfil = 'captureHiC-stratified-loops.pp.bedpe'
dis_cutoff = 80000
ratio_cutoff = 0.3
loop_set1, info_map1 = parse_loops(loop_fil1)
loop_set2, info_map2 = parse_loops(loop_fil2)
annot = []
for c in loop_set1:
    if not c in loop_set2:
        continue
    for coord in loop_set1[c]:
        loop = (c, coord[0], coord[1], c, coord[2], coord[3])

        cutoff = min(ratio_cutoff*coord[6], dis_cutoff)
        dis1 = cdist([coord[[4,5]]], loop_set2[c][:,[4,5]]).ravel()
        dis2 = cdist([coord[[0,2]]], loop_set2[c][:,[4,5]]).ravel()
        dis3 = cdist([coord[[1,3]]], loop_set2[c][:,[4,5]]).ravel()
        idx1 = np.argmin(dis1)
        idx2 = np.argmin(dis2)
        idx3 = np.argmin(dis3)
        sort_table = [
            (dis1[idx1], idx1),
            (dis2[idx2], idx2),
            (dis3[idx3], idx3)
        ]
        sort_table.sort()
        min_dis, idx = sort_table[0]
        if min_dis < cutoff:
            tmp = loop_set2[c][idx]
            key = (c, tmp[0], tmp[1], c, tmp[2], tmp[3])
            line = loop + (info_map1[loop], info_map2[key])
        else:
            line = loop + (info_map1[loop], '.')
        
        annot.append(line)
    
with open(outfil, 'w') as out:
    for tmp in annot:
        line = list(map(str, tmp))
        out.write('\t'.join(line)+'\n')