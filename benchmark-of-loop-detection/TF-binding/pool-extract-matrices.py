import numpy as np
import joblib, pyBigWig, os

bw_fils = [
    'H1ESC.EZH2.ENCFF109KCQ.bw',
    'H1ESC.POLR2A.ENCFF942TZX.bw',
    'H1ESC.CHD1.ENCFF597VKW.bw',
    'H1ESC.KDM4A.ENCFF269CHA.bw',
    'H1ESC.PHF8.ENCFF059EBB.bw',
    'H1ESC.TAF1.ENCFF837BSZ.bw',
    'H1ESC.RAD21.ENCFF056GWP.bw'
]
width = 50000
nbin_w = 50
nbin_p = 5

for bw_fil in bw_fils:
    target = bw_fil.split('.')[1]

    Indicator = '{0}.completed'.format(target)
    lockFile = '{0}.lock'.format(target)
    if os.path.exists(Indicator):
        continue
    
    if os.path.exists(lockFile):
        continue
    
    lock = open(lockFile, 'w')
    lock.close()

    bw = pyBigWig.open(bw_fil)
    labels = joblib.load('consensus-clusters.new_labels.reassigned.pkl')
    loops = joblib.load('H1ESC.umap-input.pkl')[0]
    loops_by_clusters = {}
    for i in sorted(set(labels)):
        if i == -1:
            continue
        idx = np.where(labels==i)[0]
        loops_by_clusters[i] = [loops[j][:6] for j in idx]
    
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
    
    outfil = '{0}-{1}.matrices-cache.pkl'.format('H1', target)
    joblib.dump(matrices_by_clusters, outfil, compress=('xz', 3))

    completed = open(Indicator, 'wb')
    completed.close()

    os.remove(lockFile)
