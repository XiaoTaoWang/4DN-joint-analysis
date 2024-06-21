from pyensembl import EnsemblRelease
from collections import defaultdict
import bisect, joblib
import numpy as np

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

def read_rnaseq(gene_fil, cell_type='H1'):

    genes = {}
    with open(gene_fil, 'r') as source:
        header = source.readline().rstrip().split('\t')
        cell_idx = header.index('TPM.{0}'.format(cell_type))
        for line in source:
            parse = line.rstrip().split('\t')
            gid = parse[0].split('.')[0]
            genes[gid] = float(parse[cell_idx])
    
    return genes, cell_idx

def read_imargi(gene_fil, cell_type='H1'):

    genes = {}
    with open(gene_fil, 'r') as source:
        header = source.readline().rstrip().split('\t')
        cell_idx = [i for i in range(len(header)) if header[i].startswith(cell_type)]
        for line in source:
            parse = line.rstrip().split('\t')
            gid = parse[0].split('.')[0]
            values = [float(parse[i]) for i in cell_idx]
            genes[gid] = np.mean(values)

    return genes, cell_idx

def parse_loops(loop_fil, platform='all'):

    loops = []
    with open(loop_fil, 'r') as source:
        source.readline()
        for line in source:
            c1, s1, e1, c2, s2, e2, labels = line.rstrip().split()[:7]
            labels = labels.split(',')
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            if platform != 'all':
                if not platform in labels:
                    continue
            loops.append((c1, s1, e1, c2, s2, e2))
    
    return loops

def parse_anchors(loops):

    local_anchors = defaultdict(set)
    distal_anchors = defaultdict(set)
    for line in loops:
        c1, s1, e1, c2, s2, e2 = line[:6]
        loci1 = (c1, s1, e1)
        loci2 = (c2, s2, e2)
        local_anchors[c1].add((s1, e1))
        distal_anchors[loci1].add(loci2)
        local_anchors[c2].add((s2, e2))
        distal_anchors[loci2].add(loci1)
    
    for c in local_anchors:
        local_anchors[c] = sorted(local_anchors[c])
    
    return local_anchors, distal_anchors

hg38_db = EnsemblRelease(94, species='human')
hg38_db.download()
hg38_db.index()

cell = 'HFFc6'
platform = 'all'
if cell=='H1ESC':
    cell_type = 'H1'
else:
    cell_type = 'HFF'

loops = parse_loops('{0}.union-loops.with-cluster-info.bedpe'.format(cell), platform=platform)
local_anchors, distal_anchors = parse_anchors(loops)
imargi_cache, _ = read_imargi('iMARGI_edgeR_normalized_CPM.v29.tsv', cell_type=cell_type)
rnaseq_cache, _ = read_rnaseq('RNA-seq_TPM_FPKM.v29.tsv', cell_type=cell_type)
chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']
genes = set(imargi_cache)
genes.update(set(rnaseq_cache))
anchor_by_gene = {}
for ID in genes:
    if ID in rnaseq_cache:
        tpm = rnaseq_cache[ID]
    else:
        tpm = 0
    
    if ID in imargi_cache:
        imargi = imargi_cache[ID]
    else:
        imargi = 0

    # extract distal anchors by transcripts
    gene_obj = hg38_db.gene_by_id(ID)
    c = 'chr' + gene_obj.contig
    if not c in chroms:
        continue
    s = gene_obj.start
    e = gene_obj.end
    strand = gene_obj.strand

    transcripts = gene_obj.transcripts
    tmp = {}
    for t in transcripts: 
        extracted_anchors = set()
        if strand == '+':
            key = (c, t.start-2500, t.start+2500)
        else:
            key = (c, t.end-2500, t.end+2500)
        cache = bisect_search(key, local_anchors)
        for l in cache:
            extracted_anchors.update(distal_anchors[l])
        tmp[(c, t.start, t.end)] = extracted_anchors
    
    anchor_by_gene[(c, s, e, strand, ID, tpm, imargi)] = tmp

print(len(anchor_by_gene))
joblib.dump(anchor_by_gene, '{0}.anchor-by-gene.{1}.pkl'.format(cell, platform))