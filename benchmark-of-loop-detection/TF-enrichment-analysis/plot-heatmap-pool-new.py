import numpy as np
import matplotlib, joblib, glob
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def calculate_mean(matrices):

    total = 0
    num = 0
    for M in matrices:
        total += M.sum()
        num += M.size
    
    return total/num

def sort_by_the_whole_region(cell, TFs, ID):

    cluster_ids = [1, 2, 3, 4, 5, 6]
    # determine the rank of TFs for sorting
    rank_arr = []
    m, n = 0, 1
    for TF in TFs:
        fil = '{0}-{1}.matrices-cache.pkl'.format(cell, TF)
        pool_data = joblib.load(fil)
        M1, M2 = pool_data[ID]
        value = calculate_mean([M1, M2])
        background_matrices = []
        for i in cluster_ids:
            if i != ID:
                M1, M2 = pool_data[i]
                background_matrices.extend([M1, M2])
        background_value = calculate_mean(background_matrices)
        fc = value / background_value
        rank_arr.append((fc, m, TF))
        rank_arr.append((fc, n, TF))
        m += 2
        n += 2
    
    rank_arr.sort(reverse=True)
    print([r[-1] for r in rank_arr])
    #rank_indices = [r[1] for r in rank_arr]
    rank_indices = []
    for i in range(0, len(rank_arr), 2):
        rank_indices.extend([rank_arr[i][1]-1, rank_arr[i][1]])
    print(rank_indices)

    # sort
    matrices_pool = []
    for TF in TFs:
        fil = '{0}-{1}.matrices-cache.pkl'.format(cell, TF)
        M1, M2 = joblib.load(fil)[ID]
        matrices_pool.extend([M1, M2])
    
    sort_key = []
    for i in range(M1.shape[0]):
        tmp = []
        for j in rank_indices:
            tmp.append(np.mean(matrices_pool[j][i]))
        tmp.append(i)
        sort_key.append(tuple(tmp))
    
    sort_key.sort(reverse=True)
    sort_idx = [i[-1] for i in sort_key]
    
    sort_matrices = []
    for i in range(len(matrices_pool)):
        sort_matrices.append(matrices_pool[i][sort_idx])
    
    return sort_matrices


cell = 'H1'
TFs = ['EZH2', 'POLR2A', 'CHD1', 'KDM4A', 'TAF1', 'CTCF', 'RAD21']
max_values = [18, 25, 9, 4, 12, 11, 18]
width = 50000
nbin_w = 50
nbin_p = 5
red_cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

fig = plt.figure(figsize=(7, 5))
gc = GridSpec(6, 20, figure=fig, left=0.07, right=0.95, bottom=0.05, top=0.95,
              hspace=0.1, wspace=0.4, height_ratios=[16320, 13642, 10299, 14957, 25946, 17214],
              width_ratios=[1,1,0.01,1,1,0.01,1,1,0.01,1,1,0.01,1,1,0.01,1,1,0.01,1,1])
              
col_idx = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19]
for ID in [1, 2, 3, 4, 5, 6]:
    panel_idx = [i+20*(ID-1) for i in col_idx]
    matrices_pool = sort_by_the_whole_region(cell, TFs, ID)
    for i in range(len(col_idx)//2):
        vmax = max_values[i]
        TF = TFs[i]
        # anchor1
        ti = i * 2
        M = matrices_pool[ti]
        ax1 = fig.add_subplot(gc[panel_idx[ti]])
        ax_pos1 = ax1.get_position().bounds
        sc = ax1.imshow(M, cmap=red_cmap, vmax=vmax, aspect='auto')
        if ID==6:
            ax1.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                    labelbottom=True, labeltop=False, labelleft=False, labelright=False)
            ax1.set_xticks([nbin_w+nbin_p/2])
            ax1.set_xticklabels(['anchor1'], fontsize=6)
            #ax1.set_xticks([0, nbin_w+nbin_p/2, nbin_w*2+nbin_p-1])
            #ax1.set_xticklabels(['-{0}'.format(width//1000), 'loci1', '+{0}'.format(width//1000)], fontsize=6)
        else:
            ax1.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        # anchor2
        ti = i * 2 + 1
        M = matrices_pool[ti]
        ax2 = fig.add_subplot(gc[panel_idx[ti]])
        ax_pos2 = ax2.get_position().bounds
        sc = ax2.imshow(M, cmap=red_cmap, vmax=vmax, aspect='auto')
        if ID == 6:
            ax2.tick_params(axis='both', bottom=True, top=False, left=False, right=False,
                    labelbottom=True, labeltop=False, labelleft=False, labelright=False)
            ax2.set_xticks([nbin_w+nbin_p/2])
            ax2.set_xticklabels(['anchor2'], fontsize=6)
            #ax2.set_xticks([0, nbin_w+nbin_p/2, nbin_w*2+nbin_p-1])
            #ax2.set_xticklabels(['-{0}k'.format(width//1000), 'loci2', '+{0}k'.format(width//1000)], fontsize=6)
        else:
            ax2.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        if ID==1:
            x = (ax_pos1[0] + ax_pos1[2] + ax_pos2[0]) / 2
            y = ax_pos1[1] + ax_pos1[3] + 0.015
            fig.text(x, y, TF, fontsize=7, ha='center')
        if i==0:
            x = 0.05
            y = (ax_pos1[1] + ax_pos1[1] + ax_pos1[3]) / 2
            fig.text(x, y, 'Cluster {0}'.format(ID), fontsize=7, rotation=90, ha='center', va='center')

plt.savefig('{0}-pool.svg'.format(cell), bbox_inches='tight', dpi=150)
plt.close()




