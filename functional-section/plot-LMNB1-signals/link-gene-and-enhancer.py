import joblib, bisect, pyBigWig, os, cooler, sys
from pyensembl import EnsemblRelease
from scipy.stats import gmean
import numpy as np

def parse_bed(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            if not c in D:
                D[c] = []
            D[c].append((s, e))
    
    for c in D:
        D[c].sort()
    
    return D

def bisect_search(p, ref):

    cache = set()
    if not p[0] in ref:
        return cache
    
    List = ref[p[0]]
    idx = max(0, bisect.bisect(List, p[1:])-1)
    for q in List[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add((p[0],)+q)
    
    return cache

def extract_contact_signal(M, tss, interval, res, baseline_contact=0):

    maxv = 0
    c1, s1, e1 = interval
    c2, s2, e2 = tss
    bins1 = range(s1//res, (e1+res-1)//res)
    bins2 = range(s2//res, (e2+res-1)//res)
    for b1 in bins1:
        for b2 in bins2:
            if abs(b2 - b1) > 0:
                if M[b1, b2] > maxv:
                    maxv = M[b1, b2]
    
    maxv = maxv + baseline_contact
    
    return maxv

def select_strongest_enhancer(cell, platform, res, anchor_by_gene, elements, atac, h3k27ac):

    folder = '/projects/b1100/xtwang/shares/cooler-files/4DN-joint-analysis'
    mcool = os.path.join(folder, mcool_map[cell][platform])
    clr = cooler.Cooler('{0}::resolutions/{1}'.format(mcool, res))
    if platform in ['Hi-C', 'Micro-C']:
        balance = 'weight'
    else:
        balance = False
    
    group_by_chrom = {}
    for c, s, e, strand, ID, rnaseq, imargi in anchor_by_gene:
        if not c in group_by_chrom:
            group_by_chrom[c] = []
        group_by_chrom[c].append((c, s, e, strand, ID, rnaseq, imargi))
    
    enhancer_by_ID = {}
    for c in group_by_chrom:
        M = clr.matrix(sparse=True, balance=balance).fetch(c).tocsr()
        for key in group_by_chrom[c]:
            c, s, e, strand, ID, rnaseq, imargi = key
            transcripts = anchor_by_gene[key]
            tmp_ele = set()
            for tc, ts, te in transcripts:
                distal_anchors = transcripts[(tc, ts, te)]
                for a in distal_anchors:
                    cache = bisect_search(a, elements)
                    tmp_ele.update(cache)
            
            tmp_ele = list(tmp_ele)
            if not len(tmp_ele):
                continue
                
            AC_arr = []
            for ele in tmp_ele:
                A1 = h3k27ac.stats(ele[0], ele[1], ele[2])[0]
                if A1 is None:
                    A1 = 0
                A2 = atac.stats(ele[0], ele[1], ele[2])[0]
                if A2 is None:
                    A2 = 0
                
                lens = []
                scores = []
                for tc, ts, te in transcripts:
                    lens.append(te - ts)
                    if strand == '+':
                        start = ts - 2500
                        end = ts + 2500
                    else:
                        start = te - 2500
                        end = te + 2500
                    C = extract_contact_signal(M, (tc, start, end), ele, res)
                    AC = gmean([A1, A2]) * C
                    scores.append(AC)
                
                AC_total = np.average(scores, weights=lens)
                AC_arr.append(AC_total)
            
            AC_arr = np.r_[AC_arr]
            max_i = np.argmax(AC_arr)

            enhancer_by_ID[ID] = tmp_ele[max_i]
    
    return enhancer_by_ID

mcool_map = {
    'H1ESC': {
        'Hi-C':'HIC-H1ESC-DpnII-allReps.4DNFI82R42AD.mcool',
        'Micro-C':'MicroC-H1ESC-MNase-allReps.mcool',
        'PLAC-Seq':'PLAC-H1ESC-H3K4me3-allReps.mcool',
        'ChIA-PET_CTCF':'ChIAPET-H1ESC-CTCF-allReps.mcool',
        'ChIA-PET_Pol2':'ChIAPET-H1ESC-RNAPII-allReps.mcool'
    },
    'HFFc6': {
        'Hi-C':'HIC-HFFc6-DpnII-allReps.4DNFIAVXXO55.mcool',
        'Micro-C':'MicroC-HFFc6-MNase-allReps.mcool',
        'PLAC-Seq':'PLAC-HFFc6-H3K4me3-allReps.mcool',
        'ChIA-PET_CTCF':'ChIAPET-HFFc6-CTCF-allReps.mcool',
        'ChIA-PET_Pol2':'ChIAPET-HFFc6-RNAPII-allReps.mcool'
    }
}

sig_map = {
    'H1ESC': {
        'ATAC-Seq':'H1ESC-ATAC.4DNFICPNO4M5.bw',
        'H3K27ac':'H1ESC-H3K27ac.ENCFF986PCY.bw',
        'H3K4me1':'H1ESC-H3K4me1.ENCFF088MXE.bw',
        'H3K4me3':'H1ESC-H3K4me3.ENCFF760NUN.bw',
        'H3K36me3':'H1ESC-H3K36me3.ENCFF488THD.bw'
    },
    'HFFc6': {
        'ATAC-Seq':'HFFc6-ATAC.4DNFIZ9191QU.bw',
        'H3K27ac':'HFFc6-H3K27ac.4DNFINRI6WOL.bw',
        'H3K4me1':'HFFc6-H3K4me1.4DNFI3WBAYI7.bw',
        'H3K4me3':'HFFc6-H3K4me3.4DNFIVYWZE1W.bw',
        'H3K36me3':'HFFc6-H3K36me3.ENCFF070SWD.bw'
    }
}

hg38_db = EnsemblRelease(94, species='human')
hg38_db.download()
hg38_db.index()

cell = sys.argv[1]
res = 2000
elements = parse_bed('{0}.enhancer-elements.bed'.format(cell))
anchor_by_gene = joblib.load('{0}.anchor-by-gene.all.pkl'.format(cell))
atac = pyBigWig.open(sig_map[cell]['ATAC-Seq'])
h3k27ac = pyBigWig.open(sig_map[cell]['H3K27ac'])

enhancer_by_ID = select_strongest_enhancer(cell, 'PLAC-Seq', res, anchor_by_gene, elements, atac, h3k27ac)

joblib.dump(enhancer_by_ID, '{0}.enhancer-by-gene_id.pkl'.format(cell))


