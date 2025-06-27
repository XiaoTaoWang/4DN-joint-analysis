import numpy as np
import matplotlib, joblib, os
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def generate_matrices(filenames):

    indices = []
    count = 0
    for fil in filenames:
        count += 1
        arrs_real = joblib.load(fil)
        indices.extend(list(arrs_real))

    indices_count = Counter(indices)

    indices = [i for i in indices_count if indices_count[i]==count]
    matrices_real = []
    for fil in filenames:
        arrs_real = joblib.load(fil)
        tmp = np.r_[[arrs_real[i].ravel() for i in indices]]
        matrices_real.append(tmp)

    return matrices_real

def sort_by_the_whole_region(matrices_real):

    sort_key = []
    for i in range(matrices_real[0].shape[0]):
        tmp = []
        for M in matrices_real:
            tmp.append(np.mean(M[i]))
        tmp.append(i)
        sort_key.append(tuple(tmp))
    
    sort_key.sort(reverse=True)
    sort_idx = [i[-1] for i in sort_key]

    for i in range(len(matrices_real)):
        matrices_real[i] = matrices_real[i][sort_idx]
    
    return matrices_real

my_cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
cell = 'H1ESC'

# H1ESC, 2k
filenames = [
    'H1ESC_PLAC-Seq.Enhancer_backgrounds.2kb.pkl',
    'H1ESC_ChIA-PET_Pol2.Enhancer_backgrounds.2kb.pkl',
    'H1ESC_ChIA-PET_CTCF.Enhancer_backgrounds.2kb.pkl',
    'H1ESC_MicroC.Enhancer_backgrounds.2kb.pkl',
    'H1ESC_Hi-C_FA_DSG_DpnII.Enhancer_backgrounds.2kb.pkl'
]
'''
# HFFc6, 2k
filenames = [
    'HFFc6_PLAC-Seq.Enhancer_backgrounds.2kb.pkl',
    'HFFc6_ChIA-PET_Pol2.Enhancer_backgrounds.2kb.pkl',
    'HFFc6_ChIA-PET_CTCF.Enhancer_backgrounds.2kb.pkl',
    'HFFc6_MicroC.Enhancer_backgrounds.2kb.pkl',
    'HFFc6_Hi-C_FA_DSG_DpnII.Enhancer_backgrounds.2kb.pkl'
]
'''
labels = ['H3K4me3 PLAC-Seq', 'Pol2 ChIA-PET', 'CTCF ChIA-PET',
          'Micro-C', 'Hi-C']

matrices_real = generate_matrices(filenames)
matrices_real = sort_by_the_whole_region(matrices_real)

fig = plt.figure(figsize=(6.8, 0.8))
gc = GridSpec(1, 5, figure=fig, left=0.05, right=0.95, bottom=0.05, top=0.95,
              wspace=0.3)
max_values = []
for i in range(5):
    M = matrices_real[i]
    print(labels[i], M.shape)
    ax = fig.add_subplot(gc[i])
    ax_pos = ax.get_position().bounds
    maxv = np.percentile(M.ravel(), 98)
    max_values.append(maxv)
    sc = ax.imshow(M, cmap=my_cmap, vmax=50, aspect='auto')
    ax.set_title(labels[i], fontsize=6)
    ax.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                labelbottom=True, labeltop=False, labelleft=False, labelright=False)
    ax.set_xticks([0, 15, 30])
    ax.set_xticklabels(['-30kb', 'distal\nenhancer', '+30kb'], fontsize=6)
    if i == 0:
        c_ax = fig.add_axes([ax_pos[0]-0.07, ax_pos[1]+0.7, 0.01, 0.15])
        cbar = fig.colorbar(sc, cax=c_ax)
        c_ax.tick_params(labelsize=5.5)
        cbar.set_label('distance-normalized values', size=5)

#plt.tight_layout()
plt.savefig('{0}.byenhancer.svg'.format(cell), bbox_inches='tight', dpi=150)
plt.close()