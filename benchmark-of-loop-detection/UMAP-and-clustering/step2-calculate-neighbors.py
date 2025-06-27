import os, joblib, sys
import numpy as np
from sklearnex.neighbors import NearestNeighbors

cell = sys.argv[1]
nearest_n = 3000
M, explained_variance = joblib.load('{0}.pca.pkl'.format(cell))

neigh = NearestNeighbors(n_neighbors=nearest_n,
                         algorithm='auto',
                         n_jobs=32)
neigh.fit(M)
outfil = os.path.join('{0}.{1}-nearest-neighbors.pkl'.format(cell, nearest_n))
joblib.dump(neigh, outfil, compress=('xz', 3))

