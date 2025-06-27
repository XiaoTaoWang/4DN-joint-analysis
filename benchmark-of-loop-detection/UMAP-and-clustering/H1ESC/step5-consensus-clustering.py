from collections import defaultdict, Counter
import itertools, joblib, glob, os
import igraph as ig
import numpy as np

cell = 'H1ESC'
min_node_size = 5000
strong_connection_cutoff = 0.6
res = 0.5
queue = glob.glob('community.parameter-tuning.{0}/*pkl'.format(cell))

# initialize clusters
unique_ids = defaultdict(list)
for q in queue:
    print(q)
    membership = joblib.load(q)
    for i, c in enumerate(membership):
        unique_ids[i].append(c)
    
initial_clusters = defaultdict(list)
for i in sorted(unique_ids):
    key = tuple(unique_ids[i])
    initial_clusters[key].append(i)

# build reduced graph
print('build reduced graph ...')
edges = []
weights = []
nodes = sorted(initial_clusters)
nodes_arr = np.r_[nodes]
node_ids = list(range(len(nodes)))
filtered_ids = []
remaining_ids = []
for i in node_ids:
    if len(initial_clusters[nodes[i]]) > min_node_size:
        filtered_ids.append(i)
    else:
        remaining_ids.append(i)

print('iterate within filtered clusters ...')
for i, j in itertools.combinations(filtered_ids, 2):
    mask = nodes_arr[i] == nodes_arr[j]
    ratio = mask.sum() / mask.size
    if ratio > strong_connection_cutoff: # only keep strong connections
        edges.append((i, j))
        weights.append(ratio)

print('iterate through the remaining clusters ...')
for i in remaining_ids:
    for j in filtered_ids:
        mask = nodes_arr[i] == nodes_arr[j]
        ratio = mask.sum() / mask.size
        if ratio > strong_connection_cutoff: # only keep strong connections
            edges.append((i, j))
            weights.append(ratio)

g = ig.Graph()
g.add_vertices(len(nodes))
g.add_edges(edges)
g.es['weight'] = weights
g.vs['name'] = node_ids

print('perform community detection ...')
clusters = g.community_leiden(
    weights=g.es['weight'],
    resolution=res,
    objective_function='modularity'
)
# map original nodes to final clusters
cluster_mapping = {}
for i, c in enumerate(clusters.membership):
    original = initial_clusters[nodes[i]]
    for k in original:
        cluster_mapping[k] = c
membership = []
for i in range(len(cluster_mapping)):
    membership.append(cluster_mapping[i])
counts = Counter(membership)
membership = np.r_[membership]
joblib.dump(membership, 'consensus-clusters-{0}_cutoff{1}_res{2}.pkl'.format(min_node_size, strong_connection_cutoff, res),
            compress=('xz',3))