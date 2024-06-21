import joblib, bisect, pyBigWig, os, cooler, sys
from pyensembl import EnsemblRelease
from scipy.stats import gmean
import numpy as np

def read_genes(anchor_by_gene):

    D = {}
    exp_by_gene = {}
    for key in anchor_by_gene:
        c, s, e, strand, ID, rnaseq, imargi = key
        transcripts = anchor_by_gene[key]
        for tc, ts, te in transcripts:
            if not c in D:
                D[c] = []
            if strand == '+':
                start = ts - 2500
                end = ts + 2500
            else:
                start = te - 2500
                end = te + 2500
        
            D[c].append((start, end))
            exp_by_gene[(c, start, end)] = (c, s, e, strand, ID, rnaseq, imargi)
    
    for c in D:
        D[c] =  sorted(D[c])
    
    return D, exp_by_gene

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

def extract_contact_signal(M, tss, interval, res, baseline_contact):

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

def distal_signals(cell, platform, res, anchor_by_gene, elements, atac, h3k27ac):

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

    features_by_ID = {}
    for c in group_by_chrom:
        M = clr.matrix(sparse=True, balance=balance).fetch(c).tocsr()
        d = 1000000 // res
        tmp = M.diagonal(d)
        tmp[np.isnan(tmp)] = 0
        baseline_contact = tmp.mean()
        for key in group_by_chrom[c]:
            c, s, e, strand, ID, rnaseq, imargi = key
            transcripts = anchor_by_gene[key]
            lens = []
            # Sum of interaction strength with distal elements
            C_arr = []
            # Sum of (A*C), A represents H3K27ac signals, C represents interaction strength
            AC_k27_arr = []
            # Sum of (A*C), A represents ATAC-Seq signals, C represents interaction strength
            AC_atac_arr = []
            # Sum of (A*C), A represents the geometric mean of ATAC-Seq and H3K27ac, C represents interaction strength
            AC_arr = []
            for tc, ts, te in transcripts:
                lens.append(te - ts)
                # extract proximal enhancers
                if strand == '+':
                    start = ts - 10000
                    end = ts + 10000
                else:
                    start = te - 10000
                    end = te + 10000
                
                tmp_ele = bisect_search((tc, start, end), elements)

                distal_anchors = transcripts[(tc, ts, te)]
                for a in distal_anchors:
                    cache = bisect_search(a, elements)
                    tmp_ele.update(cache)
                
                if strand == '+':
                    start = ts - 2500
                    end = ts + 2500
                else:
                    start = te - 2500
                    end = te + 2500
            
                contact_freq = 0
                AC_k27 = 0
                AC_atac = 0
                AC = 0

                for ele in tmp_ele:
                    C = extract_contact_signal(M, (tc, start, end), ele, res, baseline_contact)
                    contact_freq += C
                    A1 = h3k27ac.stats(ele[0], ele[1], ele[2])[0]
                    if A1 is None:
                        A1 = 0
                    A2 = atac.stats(ele[0], ele[1], ele[2])[0]
                    if A2 is None:
                        A2 = 0
                    
                    AC_k27 = AC_k27 + A1*C
                    AC_atac = AC_atac + A2*C
                    AC = AC + gmean([A1, A2]) * C
                
                C_arr.append(contact_freq)
                AC_k27_arr.append(AC_k27)
                AC_atac_arr.append(AC_atac)
                AC_arr.append(AC)

            arr = [
                np.average(C_arr, weights=lens),
                np.average(AC_k27_arr, weights=lens),
                np.average(AC_atac_arr, weights=lens),
                np.average(AC_arr, weights=lens)
            ]

            features_by_ID[ID] = arr
    
    return features_by_ID

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

hg38_db = EnsemblRelease(93, species='human')
hg38_db.download()
hg38_db.index()

cell = sys.argv[1]
platform = sys.argv[2]
res = 2000
elements = parse_bed('{0}.enhancer-elements.bed'.format(cell))
anchor_by_gene = joblib.load('{0}.anchor-by-gene.{1}.pkl'.format(cell, platform))
tss_by_chrom, exp_by_gene = read_genes(anchor_by_gene)
atac = pyBigWig.open(sig_map[cell]['ATAC-Seq'])
h3k27ac = pyBigWig.open(sig_map[cell]['H3K27ac'])
h3k36me3 = pyBigWig.open(sig_map[cell]['H3K36me3'])
h3k4me3 = pyBigWig.open(sig_map[cell]['H3K4me3'])
h3k4me1 = pyBigWig.open(sig_map[cell]['H3K4me1'])
pc1_hic = pyBigWig.open('{0}-PC1-25kb.Hi-C.bw'.format(cell))
pc1_microc = pyBigWig.open('{0}-PC1-25kb.Micro-C.bw'.format(cell))
pc1_h3k4me3 = pyBigWig.open('{0}-PC1-25kb.H3K4me3.bw'.format(cell))
pc1_ctcf = pyBigWig.open('{0}-PC1-25kb.CTCF.bw'.format(cell))
pc1_pol2 = pyBigWig.open('{0}-PC1-25kb.Pol2.bw'.format(cell))
tsa_rep1 = pyBigWig.open('{0}-SON_TSA-Seq.rep1.bw'.format(cell))
tsa_rep2 = pyBigWig.open('{0}-SON_TSA-Seq.rep2.bw'.format(cell))

row_info = sorted(anchor_by_gene)

# global features
print('add features of global 3d genome architectures')
global_features = []
for key in row_info:
    c, s, e, strand, ID, rnaseq, imargi = key
    lens = []
    pc1_hic_arr = []
    pc1_microc_arr = []
    pc1_h3k4me3_arr = []
    pc1_ctcf_arr = []
    pc1_pol2_arr = []
    tsa_arr_1 = []
    tsa_arr_2 = []
    transcripts = anchor_by_gene[key]
    for tc, ts, te in transcripts:
        lens.append(te - ts)
        # extract signals within [-200, 200] of a TSS
        if strand == '+':
            start = ts - 200
            end = ts + 200
        else:
            start = te - 200
            end = te + 200
    
        tmp_v = pc1_hic.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        pc1_hic_arr.append(tmp_v)

        tmp_v = pc1_microc.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        pc1_microc_arr.append(tmp_v)

        tmp_v = pc1_h3k4me3.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        pc1_h3k4me3_arr.append(tmp_v)

        tmp_v = pc1_ctcf.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        pc1_ctcf_arr.append(tmp_v)

        tmp_v = pc1_pol2.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        pc1_pol2_arr.append(tmp_v)

        tmp_v = tsa_rep1.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        tsa_arr_1.append(tmp_v)

        tmp_v = tsa_rep2.stats(c, start, end)[0]
        if tmp_v is None:
            tmp_v = 0
        tsa_arr_2.append(tmp_v)

    pc1_hic_score = np.average(pc1_hic_arr, weights=lens)
    pc1_microc_score = np.average(pc1_microc_arr, weights=lens)
    pc1_h3k4me3_score = np.average(pc1_h3k4me3_arr, weights=lens)
    pc1_ctcf_score = np.average(pc1_ctcf_arr, weights=lens)
    pc1_pol2_score = np.average(pc1_pol2_arr, weights=lens)
    tsa_score_1 = np.average(tsa_arr_1, weights=lens)
    tsa_score_2 = np.average(tsa_arr_2, weights=lens)
    global_features.append([pc1_hic_score, pc1_microc_score, pc1_h3k4me3_score, pc1_ctcf_score, pc1_pol2_score, tsa_score_1, tsa_score_2])

y = []
# ChIP/ATAC-seq signals around TSS
print('ChIP/ATAC-seq signals around TSS')
local_features = []
for key in row_info:
    c, s, e, strand, ID, rnaseq, imargi = key
    y.append([rnaseq, imargi])
    lens = []
    k4me3_arr = []
    atac_arr = []
    ac_arr = []
    k4me1_arr = []
    transcripts = anchor_by_gene[key]
    for tc, ts, te in transcripts:
        lens.append(te - ts)
        # extract signals within [-200, 200] of a TSS
        if strand == '+':
            start = ts - 200
            end = ts + 200
        else:
            start = te - 200
            end = te + 200
    
        local_k4me3 = h3k4me3.stats(c, start, end)[0]
        if local_k4me3 is None:
            local_k4me3 = 0
        k4me3_arr.append(local_k4me3)
        local_atac = atac.stats(c, start, end)[0]
        if local_atac is None:
            local_atac = 0
        atac_arr.append(local_atac)
        local_ac = h3k27ac.stats(c, start, end)[0]
        if local_ac is None:
            local_ac = 0
        ac_arr.append(local_ac)
        local_k4me1 = h3k4me1.stats(c, start, end)[0]
        if local_k4me1 is None:
            local_k4me1 = 0
        k4me1_arr.append(local_k4me1)

    local_k4me3 = np.average(k4me3_arr, weights=lens)
    local_atac = np.average(atac_arr, weights=lens)
    local_ac = np.average(ac_arr, weights=lens)
    local_k4me1 = np.average(k4me1_arr, weights=lens)
    local_k36 = h3k36me3.stats(c, s, e)[0] # whole gene body
    if local_k36 is None:
        local_k36 = 0
    local_features.append([local_k4me3, local_atac, local_ac, local_k4me1, local_k36])

# number of distal elements
print('number of distal enhancers')
num_features = []
for key in row_info:
    c, s, e, strand, ID, rnaseq, imargi = key
    transcripts = anchor_by_gene[key]
    lens = []
    nums = []
    for tc, ts, te in transcripts:
        lens.append(te - ts)
        # extract proximal elements
        if strand == '+':
            start = ts - 10000
            end = ts + 10000
        else:
            start = te - 10000
            end = te + 10000
        
        tmp_ele = bisect_search((tc, start, end), elements)
        distal_anchors = transcripts[(tc, ts, te)]
        for a in distal_anchors:
            cache = bisect_search(a, elements)
            tmp_ele.update(cache)
        nums.append(len(tmp_ele))
    
    num_features.append([np.average(nums, weights=lens)])

# expression of distal genes
exp_features = []
for key in row_info:
    transcripts = anchor_by_gene[key]
    lens = []
    exp_max_arr = []
    exp_avg_arr = []
    for tc, ts, te in transcripts:
        lens.append(te - ts)
        distal_anchors = transcripts[(tc, ts, te)]
        distal_genes = set()
        for a in distal_anchors:
            cache = bisect_search(a, tss_by_chrom)
            for t in cache:
                distal_genes.add(exp_by_gene[t])

        if len(distal_genes):
            rnaseq_pool = []
            gene_lens = []
            for c, s, e, strand, ID, rnaseq, imargi in distal_genes:
                rnaseq_pool.append(rnaseq)
                gene_lens.append(e - s)
            exp_max_arr.append(max(rnaseq_pool))
            exp_avg_arr.append(np.average(rnaseq_pool, weights=gene_lens))
        else:
            exp_max_arr.append(0)
            exp_avg_arr.append(0)
    
    exp_features.append([np.average(exp_max_arr, weights=lens), np.average(exp_avg_arr, weights=lens)])

# ATAC-Seq/ChIP-Seq signals at distal elements
distal_activities = []
for key in row_info:
    c, s, e, strand, ID, rnaseq, imargi = key
    transcripts = anchor_by_gene[key]
    lens = []
    score_arr1 = []
    score_arr2 = []
    score_arr3 = []
    for tc, ts, te in transcripts:
        lens.append(te - ts)
        # extract proximal elements
        if strand == '+':
            start = ts - 10000
            end = ts + 10000
        else:
            start = te - 10000
            end = te + 10000
        
        atac_score = 0
        ac_score = 0
        atac_ac_score = 0
        tmp_ele = bisect_search((tc, start, end), elements)

        distal_anchors = transcripts[(tc, ts, te)]
        for a in distal_anchors:
            cache = bisect_search(a, elements)
            tmp_ele.update(cache)

        for c, s, e in tmp_ele:
            e_atac = atac.stats(c, s, e)[0]
            if e_atac is None:
                e_atac = 0
            e_ac = h3k27ac.stats(c, s, e)[0]
            if e_ac is None:
                e_ac = 0
            
            atac_score += e_atac
            ac_score += e_ac
            atac_ac_score += gmean([e_atac, e_ac])
        
        score_arr1.append(atac_score)
        score_arr2.append(ac_score)
        score_arr3.append(atac_ac_score)
    
    score1 = np.average(score_arr1, weights=lens)
    score2 = np.average(score_arr2, weights=lens)
    score3 = np.average(score_arr3, weights=lens)
    distal_activities.append([score1, score2, score3])

# consider chromatin contacts
print('Hi-C')
# Hi-C
hic = distal_signals(cell, 'Hi-C', res, anchor_by_gene, elements, atac, h3k27ac)
hic_features = []
for key in row_info:
    hic_features.append(hic[key[4]])

# Micro-C
print('Micro-C')
microc = distal_signals(cell, 'Micro-C', res, anchor_by_gene, elements, atac, h3k27ac)
microc_features = []
for key in row_info:
    microc_features.append(microc[key[4]])

# PLAC-Seq
print('PLAC-Seq')
plac = distal_signals(cell, 'PLAC-Seq', res, anchor_by_gene, elements, atac, h3k27ac)
plac_features = []
for key in row_info:
    plac_features.append(plac[key[4]])

# CTCF
print('CTCF')
ctcf = distal_signals(cell, 'ChIA-PET_CTCF', res, anchor_by_gene, elements, atac, h3k27ac)
ctcf_features = []
for key in row_info:
    ctcf_features.append(ctcf[key[4]])

# pol2
print('Pol2')
pol2 = distal_signals(cell, 'ChIA-PET_Pol2', res, anchor_by_gene, elements, atac, h3k27ac)
pol2_features = []
for key in row_info:
    pol2_features.append(pol2[key[4]])

print('output matrix ...')
fea_matrix = []
for i in range(len(pol2_features)):
    fea_matrix.append(local_features[i] + num_features[i] + exp_features[i] + distal_activities[i] + hic_features[i] + microc_features[i] + plac_features[i] + ctcf_features[i] + pol2_features[i] + global_features[i])

fea_matrix = np.r_[fea_matrix]
y = np.r_[y]

feature_labels = [
    'tss_h3k4me3',
    'tss_atac',
    'tss_h3k27ac',
    'tss_h3k4me1',
    'local_h3k36me3',
    'num_distal_enhancers',
    'distal_rnaseq_max',
    'distal_rnaseq_avg',
    'A_atac',
    'A_h3k27ac',
    'A_atac_h3k27ac'
]
platform_specific = [
    '{0}_C',
    '{0}_A*C_h3k27ac',
    '{0}_A*C_atac',
    '{0}_A*C_atac_h3k27ac'
]
for p in ['Hi-C', 'Micro-C', 'H3K4me3', 'CTCF', 'Pol2']:
    for pattern in platform_specific:
        feature_labels.append(pattern.format(p))
feature_labels.extend(['Hi-C_compartment',
                       'Micro-C_compartment',
                       'H3K4me3_compartment',
                       'CTCF_compartment',
                       'Pol2_compartment',
                       'SON_TSA-Seq_rep1',
                       'SON_TSA-Seq_rep2'])

joblib.dump([fea_matrix, y, row_info, feature_labels],
            '{0}-prediction-features.{1}.{2}k.pkl'.format(cell, platform, res//1000),
            compress=('xz', 3))

