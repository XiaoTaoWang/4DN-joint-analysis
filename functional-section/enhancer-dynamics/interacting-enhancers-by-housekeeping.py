import sys, os, bisect, joblib, matplotlib
import numpy as np
from collections import defaultdict, Counter
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import OrRd_4, OrRd_9
from palettable.colorbrewer.diverging import RdYlGn_4_r
import matplotlib.gridspec as gridspec
from palettable.colorbrewer.sequential import Greys_3
import seaborn as sns

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

hg38_db = EnsemblRelease(94, species='human')
hg38_db.download()
hg38_db.index()
cmap = RdYlGn_4_r.mpl_colormap

def parse_loops(loop_fil):

    loops = []
    loops_by_anchors = {}
    with open(loop_fil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2, labels = line.rstrip().split()
            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
            loops.append((c1, s1, e1, c2, s2, e2))
            key1 = (c1, s1, e1)
            key2 = (c2, s2, e2)
            if not key1 in loops_by_anchors:
                loops_by_anchors[key1] = []
            if not key2 in loops_by_anchors:
                loops_by_anchors[key2] = []
            loops_by_anchors[key1].append((c1, s1, e1, c2, s2, e2, labels))
            loops_by_anchors[key2].append((c1, s1, e1, c2, s2, e2, labels))
    
    return loops, loops_by_anchors

def parse_enhancers(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom, start, end = parse[:3]
            start, end = int(start), int(end)
            if not chrom in D:
                D[chrom] = set()
            D[chrom].add((start, end))
    
    for c in D:
        D[c] = sorted(D[c])
    
    return D

def read_genes(gene_fil, hg38_db):

    tss = []
    chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
    with open(gene_fil, 'r') as source:
        source.readline()
        for line in source:
            gname = line.rstrip().split(';')[1]
            try:
                tmp = hg38_db.genes_by_name(gname)
                if len(tmp) == 1:
                    ID = tmp[0].gene_id
                    chrom, start, end, strand = 'chr'+tmp[0].contig, tmp[0].start, tmp[0].end, tmp[0].strand
                    if not chrom in chroms:
                        continue
                    if strand == '+':
                        tss.append((chrom, start-2500, start+2500, strand, ID))
                    else:
                        tss.append((chrom, end-2500, end+2500, strand, ID))
            except:
                pass
    
    return tss

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

def parse_and_filter_anchors(loops, enhancers):

    distal_anchors = defaultdict(set)
    distal_enhancers = defaultdict(set)
    for line in loops:
        c1, s1, e1, c2, s2, e2 = line[:6]
        loci1 = (c1, s1, e1)
        loci2 = (c2, s2, e2)
        cache1 = bisect_search(loci1, enhancers)
        cache2 = bisect_search(loci2, enhancers)
        if len(cache1):
            distal_anchors[c2].add((s2, e2))
            for ec, es, ee in cache1:
                distal_enhancers[loci2].add((ec, es, ee))
        if len(cache2):
            distal_anchors[c1].add((s1, e1))
            for ec, es, ee in cache2:
                distal_enhancers[loci1].add((ec, es, ee))
    
    for c in distal_anchors:
        distal_anchors[c] = sorted(distal_anchors[c])
    
    return distal_anchors, distal_enhancers

cell = 'H1ESC'

loops, loops_by_anchors = parse_loops('../{0}.union-loops.bedpe'.format(cell))

enhancers = parse_enhancers('{0}.enhancer-elements.bed'.format(cell))
distal_anchors, distal_enhancers = parse_and_filter_anchors(loops, enhancers)
tss = read_genes('Housekeeping_GenesHuman.csv', hg38_db)

print('total house-keeping genes: {0}'.format(len(tss)))

enhancer_count = defaultdict(list)
house_keeping_classes = {}
extracted_loops = []
house_keeping_list = []
for c, s, e, strand, ID in tss:
    extracted_anchors = set()
    gene_obj = hg38_db.gene_by_id(ID)
    transcripts = gene_obj.transcripts
    for t in transcripts:
        if strand == '+':
            key = (c, t.start-2500, t.start+2500)
        else:
            key = (c, t.end-2500, t.end+2500)
        cache = bisect_search(key, distal_anchors)
        extracted_anchors.update(cache)
    extracted_enhancers = set()
    for loci in extracted_anchors:
        extracted_enhancers.update(distal_enhancers[loci])
        extracted_loops.extend(loops_by_anchors[loci])
    extracted_enhancers = list(extracted_enhancers)

    house_keeping_classes[ID] = len(extracted_enhancers)
    house_keeping_list.append([ID, gene_obj.gene_name])

with open('{0}.human-housekeeping.txt'.format(cell), 'w') as out:
    for ID, name in house_keeping_list:
        out.write('\t'.join([ID, name, str(house_keeping_classes[ID])])+'\n')

## house-keeping genes
colors = OrRd_9.hex_colors
collect = list(house_keeping_classes.values())
collect_counts = Counter(collect)
total = len(house_keeping_classes)
ratios = [
    collect_counts[0] / total,
    collect_counts[1] / total,
    collect_counts[2] / total,
    collect_counts[3] / total,
    collect_counts[4] / total,
    collect_counts[5] / total,
    sum([collect_counts[i] for i in collect_counts if 5 < i < 8]) / total,
    sum([collect_counts[i] for i in collect_counts if 7 < i < 11]) / total,
    sum([collect_counts[i] for i in collect_counts if i > 10]) / total
]
fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot(111)
patches, texts, autotexts = ax.pie(ratios[::-1], colors=colors[::-1], autopct='%1.1f%%',
                                startangle=90, labeldistance=1.03, pctdistance=0.6)
plt.savefig('{0}.test.house-keeping.enhancer-counts.svg'.format(cell), dpi=300, bbox_inches='tight')
plt.close()
