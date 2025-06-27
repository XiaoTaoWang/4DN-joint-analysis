import joblib, os, sys
import numpy as np
from igraph import Graph

def flip_features(ori_loops, ori_fea):

    loops = []
    fea = []
    for i in range(len(ori_loops)):
        c1, s1, e1, c2, s2, e2, label = ori_loops[i]
        arr = ori_fea[i]
        anchor1 = arr[:11]
        anchor2 = arr[11:]
        if anchor1.max() < anchor2.max():
            anchor1, anchor2 = anchor2, anchor1

        loops.append((c1, s1, e1, c2, s2, e2, label))
        fea.append(np.r_[anchor1, anchor2])
    
    fea = np.r_[fea]

    return loops, fea

cell = sys.argv[1]
infil = '{0}.3000-nearest-neighbors.pkl'.format(cell)
folder = 'community.parameter-tuning.{0}'.format(cell)
if not os.path.exists(folder):
    os.mkdir(folder)

ks = [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90] + list(range(100, 3000, 100))
resolution_values = [0.5, 1, 1.5]
for k in ks:
    for res in resolution_values:
        Indicator = '{0}-{1}-{2:.4g}.completed'.format(cell, k, res)
        lockFile = '{0}-{1}-{2:.4g}.lock'.format(cell, k, res)
        if os.path.exists(Indicator):
            continue

        if os.path.exists(lockFile):
            continue

        lock = open(lockFile, 'w')
        lock.close()

        neigh = joblib.load(infil)
        M = neigh.kneighbors_graph(n_neighbors=k, mode='connectivity')
        x, y = M.nonzero()
        g = Graph(n=M.shape[0], edges=list(zip(x, y)), directed=False)
        clusters = g.community_leiden( 
            resolution=res,
            objective_function='modularity',
        )
        membership = np.r_[clusters.membership]
        outfil = os.path.join(folder, 'leiden-membership_k{0}_res{1:.4g}.pkl'.format(k, res))
        joblib.dump(clusters.membership, outfil, compress=('xz', 3))

        completed = open(Indicator, 'w')
        completed.close()

        os.remove(lockFile)

