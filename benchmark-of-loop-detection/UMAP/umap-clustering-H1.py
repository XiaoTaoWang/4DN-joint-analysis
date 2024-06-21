import numpy as np
import matplotlib.pyplot as plt
import umap, joblib
import umap.plot
import hdbscan
from collections import Counter

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

h1_loops, h1_fea = joblib.load('H1ESC.umap-input.pkl')
h1_loops, h1_fea = flip_features(h1_loops, h1_fea)
hff_loops, hff_fea = joblib.load('HFFc6.umap-input.pkl')
hff_loops, hff_fea = flip_features(hff_loops, hff_fea)

loops = h1_loops + hff_loops
fea = np.r_[h1_fea, hff_fea]
print(fea.shape)
'''
# step 1
# parameter test
for n_n in [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]:
    mapper = umap.UMAP(n_neighbors=n_n, min_dist=0, n_components=2, metric='correlation', random_state=42).fit(fea)
    joblib.dump(mapper, '{0}.correlation-{1}.pkl'.format('union', n_n), compress=('xz', 3))
    for m_s in [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]:
        for m_c in [400, 500, 600]:
            print(n_n, m_s, m_c)
            labels = hdbscan.HDBSCAN(min_samples=m_s, min_cluster_size=500).fit_predict(mapper.embedding_)
            umap.plot.points(mapper, labels=labels, show_legend=True)
            plt.savefig('{0}.correlation-{1}.hdbscan{2}_{3}.png'.format('union', n_n, m_s, m_c), dpi=1000)
            plt.close()
# step 2, fine tune the HDBSCAN parameters
for n_n in [40, 50]:
    mapper = joblib.load('{0}.correlation-{1}.pkl'.format('union', n_n))
    for m_s in range(30, 50):
        for m_c in range(300, 1050, 50):
            print(n_n, m_s, m_c)
            labels = hdbscan.HDBSCAN(min_samples=m_s, min_cluster_size=500).fit_predict(mapper.embedding_)
            umap.plot.points(mapper, labels=labels, show_legend=True)
            plt.savefig('{0}.correlation-{1}.hdbscan{2}_{3}.png'.format('test', n_n, m_s, m_c), dpi=1000)
            plt.close()

# step 3
mapper = joblib.load('union.correlation-40.pkl')
labels = hdbscan.HDBSCAN(min_samples=47, min_cluster_size=1000).fit_predict(mapper.embedding_)
color_key = {
    -1:'#e1e1e1', 0:'#e1e1e1', 1:'#e1e1e1', 3:'#e1e1e1',
    2: '#1F78B4', 4:'#E31A1C', 5:'#B2DF8A', 6:'#A6CEE3',
    7:'#FB9A99', 8:'#FDBF6F'
}
ax = umap.plot.points(mapper, labels=labels, show_legend=True, color_key=color_key)
ax.set_xlim(260, 600)
ax.set_ylim(230, 550)
plt.savefig('test.svg', dpi=1000)
plt.close()

# step 3
cell_labels = np.r_[np.ones(len(h1_loops), dtype=int),
                    np.zeros(len(hff_loops), dtype=int)]
ax = umap.plot.points(mapper, labels=cell_labels, color_key={0:'#e0e0e0', 1:'#E31A1C'}, show_legend=True)
ax.set_xlim(260, 600)
ax.set_ylim(230, 550)
plt.savefig('corr40-6clusters.H1.svg', dpi=1000)
plt.close()

cell_labels = np.r_[np.zeros(len(h1_loops), dtype=int),
                    np.ones(len(hff_loops), dtype=int)]
ax = umap.plot.points(mapper, labels=cell_labels, color_key={0:'#e0e0e0', 1:'#E31A1C'}, show_legend=True)
ax.set_xlim(260, 600)
ax.set_ylim(230, 550)
plt.savefig('corr40-6clusters.HFF.svg', dpi=1000)
plt.close()
'''
# step 4
mapper = joblib.load('union.correlation-40.pkl')
labels = hdbscan.HDBSCAN(min_samples=47, min_cluster_size=1000).fit_predict(mapper.embedding_)
color_key = {
    -1:'#ffffff', 0:'#ffffff', 1:'#ffffff', 3:'#ffffff',
    2: '#1F78B4', 4:'#E31A1C', 5:'#B2DF8A', 6:'#A6CEE3',
    7:'#FB9A99', 8:'#FDBF6F'
}
h1_labels = labels.copy()
h1_labels[len(h1_loops):] = -1
ax = umap.plot.points(mapper, labels=h1_labels, show_legend=True, color_key=color_key)
ax.set_xlim(260, 600)
ax.set_ylim(230, 550)
plt.savefig('H1-umap-6clusters.svg', dpi=1000)
plt.close()

hff_labels = labels.copy()
hff_labels[:len(h1_loops)] = -1
ax = umap.plot.points(mapper, labels=hff_labels, show_legend=True, color_key=color_key)
ax.set_xlim(260, 600)
ax.set_ylim(230, 550)
plt.savefig('HFF-umap-6clusters.svg', dpi=1000)
plt.close()

platforms = ['Hi-C', 'Micro-C', 'PLACSeq-H3K4me3', 'ChIAPET-CTCF', 'ChIAPET-RNAPII']
for t in platforms:
    platform_labels = []
    for l in h1_loops:
        if t in l[-1]:
            platform_labels.append(1)
        else:
            platform_labels.append(0)
    platform_labels = np.r_[platform_labels, np.zeros(len(hff_loops), dtype=int)]
    ax = umap.plot.points(mapper, labels=platform_labels, color_key={0:'#ffffff', 1:'#E31A1C'}, show_legend=True)
    ax.set_xlim(260, 600)
    ax.set_ylim(230, 550)
    plt.savefig('corr40-6clusters.H1.{0}.svg'.format(t), dpi=1000)
    plt.close()

    platform_labels = []
    for l in hff_loops:
        if t in l[-1]:
            platform_labels.append(1)
        else:
            platform_labels.append(0)
    platform_labels = np.r_[np.zeros(len(h1_loops), dtype=int), platform_labels]
    ax = umap.plot.points(mapper, labels=platform_labels, color_key={0:'#ffffff', 1:'#E31A1C'}, show_legend=True)
    ax.set_xlim(260, 600)
    ax.set_ylim(230, 550)
    plt.savefig('corr40-6clusters.HFF.{0}.svg'.format(t), dpi=1000)
    plt.close()

