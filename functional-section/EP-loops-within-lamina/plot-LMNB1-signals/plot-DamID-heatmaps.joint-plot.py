import joblib, bisect, matplotlib, pyBigWig
import numpy as np
from collections import defaultdict, Counter
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from pyensembl import EnsemblRelease
from matplotlib.colors import Normalize

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

hg38_db = EnsemblRelease(94, species='human')
hg38_db.download()
hg38_db.index()

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        if self.vmin == self.vmax:
            return np.ma.masked_array(np.interp(value, [self.vmin], [0.5]))
        
        if self.vmin < self.midpoint < self.vmax:
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        elif self.vmin >= self.midpoint:
            x, y = [self.vmin, self.vmax], [0.5, 1]
        elif self.vmax <= self.midpoint:
            x, y = [self.vmin, self.vmax], [0, 0.5]
            
        return np.ma.masked_array(np.interp(value, x, y))

def extract_bigwig_values(gene_list, damid_fil, enhancers_by_ID, width=100000,
    n_gene_bins=10, n_neighbor_bins=20):

    L = []
    sort_key = []
    db = pyBigWig.open(damid_fil)
    chromosizes = db.chroms()
    idx = 0
    for g in gene_list:
        c, s, e, strand, gene_id, tpm = g
        gene_length = e - s
        if not enhancers_by_ID is None:
            if not gene_id in enhancers_by_ID:
                continue

        # gene    
        left = s - width
        right = e + width
        if left < 0:
            continue
        if right > chromosizes[c]:
            continue

        arr1 = db.stats(c, s, e, nBins=n_gene_bins)
        for i in range(len(arr1)):
            if arr1[i] is None:
                arr1[i] = 0

        arr2 = db.stats(c, left, s, nBins=n_neighbor_bins)
        for i in range(len(arr2)):
            if arr2[i] is None:
                arr2[i] = 0

        arr3 = db.stats(c, e, right, nBins=n_neighbor_bins)
        for i in range(len(arr3)):
            if arr3[i] is None:
                arr3[i] = 0

        arr_gene = arr2 + arr1 + arr3
        if strand=='-':
            arr_gene = arr_gene[::-1]
        arr_gene = np.r_[arr_gene]

        # enhancer
        enhancer_loci = enhancers_by_ID[gene_id]
        left = enhancer_loci[1] - width
        right = enhancer_loci[2] + width
        if left < 0:
            continue
        if right > chromosizes[c]:
            continue

        arr_enhancer = db.stats(c, left, right, nBins=51)
        for i in range(len(arr_enhancer)):
            if arr_enhancer[i] is None:
                arr_enhancer[i] = 0
        arr_enhancer = np.r_[arr_enhancer]

        L.append((arr_gene, arr_enhancer, tpm, gene_length))
        if strand=='+':
            sort_key.append((arr1[0], idx))
        else:
            sort_key.append((arr1[-1], idx))
        idx += 1
    
    sort_key.sort()
    idx_queue = [i[1] for i in sort_key]
    M_gene = np.r_[[L[i][0] for i in idx_queue]]
    M_enhancer = np.r_[[L[i][1] for i in idx_queue]]
    exp_arr = np.r_[[L[i][2] for i in idx_queue]]
    gene_lens = np.r_[[L[i][3] for i in idx_queue]]

    return M_gene, M_enhancer, exp_arr, gene_lens

cell = 'HFFc6'
gene_fil = 'HFFc6.lamina-genes.pkl'
damid_fil = 'HFFc6-DamID.4DNFI7724Y7Q.bw'
'''
cell = 'H1ESC'
gene_fil = 'H1ESC.lamina-genes.pkl'
damid_fil = 'H1ESC-DamID.4DNFI6BH48Y3.bw'
'''
enhancers_by_ID = joblib.load('{0}.enhancer-by-gene_id.pkl'.format(cell))
_, genes_with_enhancers = joblib.load(gene_fil)
M_gene, M_enhancer, exp, gene_lens  = extract_bigwig_values(genes_with_enhancers, damid_fil, enhancers_by_ID)

fig = plt.figure(figsize=(3.5, 4))
gc = GridSpec(2, 4, figure=fig, left=0.1, right=0.9, bottom=0.1, top=0.9,
              hspace=0.03, wspace=0.3, #height_ratios=[M_enhancer.shape[0], M_no.shape[0]],
              height_ratios=[3,1], width_ratios=[3.5, 3.5, 1, 1])

# normalization core for DamID
vmin = min(M_gene.min(), M_enhancer.min())
vmax = max(M_gene.max(), M_enhancer.max())
print(vmin, vmax)
norm = MidpointNormalize(vmin=-1, vmax=2.5, midpoint=0)

# gene DamID-Seq
ax1 = fig.add_subplot(gc[0])
sc = ax1.imshow(M_gene, cmap='RdBu_r', norm=norm, aspect='auto')
ax1.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                labelbottom=True, labeltop=False, labelleft=False, labelright=False)
ax1.set_xticks([0, 20, 30, 50])
ax1.set_xticklabels(['-100kb', 'TSS', 'TES', '+100kb'], fontsize=6)

ax_pos = ax1.get_position().bounds
x = 0.03
c_ax = fig.add_axes([x, ax_pos[1]+0.5, 0.03, 0.1])
cbar = fig.colorbar(sc, cax=c_ax, ticks=[-1, 0, 2.5])
c_ax.tick_params(labelsize=6)
cbar.set_label('LMNB1 DamID-Seq\nnormalized counts', size=6)

# enhancer DamID-Seq
ax2 = fig.add_subplot(gc[1])
sc = ax2.imshow(M_enhancer, cmap='RdBu_r', norm=norm, aspect='auto')
ax2.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                labelbottom=True, labeltop=False, labelleft=False, labelright=False)
ax2.set_xticks([0, 25, 50])
ax2.set_xticklabels(['-100kb', 'distal\nenhancer', '+100kb'], fontsize=6)

# RNA-Seq
ax3 = fig.add_subplot(gc[2])
ax3.fill_betweenx(range(len(exp)), x1=0, x2=exp[::-1], color='#E41A1C', edgecolor='none', alpha=0.8)
ax3.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                labelbottom=True, labeltop=False, labelleft=False, labelright=False)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_xlabel('RNA log2(TPM+1)', fontsize=6)
ax3.set_ylim(0, len(exp)-1)
ax3.set_xlim(0, 5)

# gene lengths
ax4 = fig.add_subplot(gc[3])
ax4.fill_betweenx(range(len(gene_lens)), x1=0, x2=gene_lens/1000, color='#999999', edgecolor='none', alpha=0.9)
ax4.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                labelbottom=True, labeltop=False, labelleft=False, labelright=False)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.set_xlabel('Gene Length (kb)', fontsize=6)
ax4.set_ylim(0, len(gene_lens)-1)
#ax4.set_xlim(0, 5)
    
for ax in [ax1, ax2, ax3, ax4]:
    ax.xaxis.set_tick_params(width=1, labelsize=6, pad=2)
    ax.yaxis.set_tick_params(width=1, labelsize=6, pad=2)

plt.savefig('{0}-DamID-heatmaps.joint.svg'.format(cell), bbox_inches='tight', dpi=150)
plt.close()