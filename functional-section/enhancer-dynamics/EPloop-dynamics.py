import joblib, bisect
from collections import Counter, defaultdict
import matplotlib
import matplotlib.pyplot as plt
from pybedtools import BedTool

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

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

def overlap(a1, a2):

    t1 = list(a1[1:])
    t2 = list(a2[1:])
    if (t1[1] <= t2[0]) or (t2[1] <= t1[0]):
        return 0
    
    mi = t1 + t2
    mi.sort()
    OR = (mi[2]-mi[1])/(mi[3]-mi[0]) # intersect / union

    return OR

def calculate_union_intervals(h1_anchors, hff_anchors):

    intervals = []
    for a in h1_anchors:
        intervals.append(a)
    
    for a in hff_anchors:
        intervals.append(a)
    
    intervals = sorted(intervals)
    bed = BedTool(intervals)
    merged = bed.merge()
    union = []
    for r in merged:
        union.append((r.chrom, r.start, r.end))

    return union

h1_enhancers = parse_bed('H1ESC.enhancer-elements.bed')
hff_enhancers = parse_bed('HFFc6.enhancer-elements.bed')
h1_anchor_by_gene = joblib.load('H1ESC.housekeeping.interacting-loci.pkl')
hff_anchor_by_gene = joblib.load('HFFc6.housekeeping.interacting-loci.pkl')

EP_loop_info = [] # binary vectors, columns: H1 loop, HFF loop, H1 enhancer, and HFF enhancer
loop_intermediate = defaultdict(set)

for ID in h1_anchor_by_gene:
    h1_anchors = h1_anchor_by_gene[ID]
    hff_anchors = hff_anchor_by_gene[ID]
    union = calculate_union_intervals(h1_anchors, hff_anchors)
    annot = {}
    for u_a in union:
        h1_score1 = 0
        hff_score1 = 0
        for h1_a in h1_anchors:
            OR = overlap(u_a, h1_a)
            if OR > h1_score1:
                h1_score1 = OR
        
        for hff_a in hff_anchors:
            OR = overlap(u_a, hff_a)
            if OR > hff_score1:
                hff_score1 = OR
        
        if (h1_score1 > 0) and (hff_score1 > 0):
            annot[u_a] = [1, 1]
        elif (h1_score1 > 0) and (hff_score1 == 0):
            annot[u_a] = [1, 0]
        elif (h1_score1 == 0) and (hff_score1 > 0):
            annot[u_a] = [0, 1]
        else:
            annot[u_a] = [0, 0]
    
    # whether there are enhancers within that anchor in H1/HFF
    for a in annot:
        h1_cache = bisect_search(a, h1_enhancers)
        hff_cache = bisect_search(a, hff_enhancers)
        if len(h1_cache) and len(hff_cache):
            annot[a].extend([1, 1])
        elif len(h1_cache) and (not len(hff_cache)):
            annot[a].extend([1, 0])
        elif (not len(h1_cache)) and len(hff_cache):
            annot[a].extend([0, 1])
        else:
            annot[a].extend([0, 0])
        
        key = tuple(annot[a])
        c, ts, te = ID[:3]
        c, es, ee = a
        if es > te:
            loop_intermediate[key].add((c, ts, te, c, es, ee))
        else:
            loop_intermediate[key].add((c, es, ee, c, ts, te))
    
    for a in annot:
        EP_loop_info.append(tuple(annot[a]))

counts = Counter(EP_loop_info)
# pie chart
count1 = counts[(1, 1, 1, 1)] # shared loops, shared enhancers
count2 = counts[(1, 1, 1, 0)] + counts[(1, 1, 0, 1)] # shared loops, different enhancers
count3 = counts[(1, 0, 1, 1)] + counts[(0, 1, 1, 1)] # different loops, shared enhancers
count4 = counts[(1, 0, 1, 0)] + counts[(0, 1, 0, 1)] # different loops, different enhancers
colors = ['#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4']
total = count1 + count2 + count3 + count4
ratios = [
    count1 / total,
    count2 / total,
    count3 / total,
    count4 / total
]
fig = plt.figure(figsize=(2.5, 2.5))
ax = fig.add_subplot(111)
patches, texts, autotexts = ax.pie(ratios[::-1], colors=colors[::-1], autopct='%1.1f%%',
                                startangle=90, labeldistance=1.03, pctdistance=0.6)
plt.savefig('EP-loop-dynamics.piechart.svg', dpi=300, bbox_inches='tight')
plt.close()

name_map = {
    (1, 1, 1, 1): 'both_shared.bedpe',
    (1, 1, 1, 0): 'shared_loops.h1_specific_enhancers.bedpe',
    (1, 1, 0, 1): 'shared_loops.hff_specific_enhancers.bedpe',
    (1, 0, 1, 1): 'h1_specific_loops.shared_enhancers.bedpe',
    (0, 1, 1, 1): 'hff_specific_loops.shared_enhancers.bedpe',
    (1, 0, 1, 0): 'h1_specific_loops.h1_specific_enhancers.bedpe',
    (0, 1, 0, 1): 'hff_specific_loops.hff_specific_enhancers.bedpe'
}
for k in name_map:
    with open(name_map[k], 'w') as out:
        for loop in sorted(loop_intermediate[k]):
            out.write('\t'.join(list(map(str, loop)))+'\n')