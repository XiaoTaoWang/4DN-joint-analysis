import joblib, sys, os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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
loops, fea = joblib.load('{0}.umap-input.pkl'.format(cell))
loops, fea = flip_features(loops, fea)

pca = PCA()
pca.fit(fea)
explained_variance = pca.explained_variance_ratio_
print(explained_variance)
fea_new = pca.transform(fea)
joblib.dump([fea_new, explained_variance], '{0}.pca.pkl'.format(cell), compress=('xz', 3))

