import os, joblib, glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

new_rc_params = {
    'text.usetex': False,
    'svg.fonttype': 'none'
}

matplotlib.rcParams.update(new_rc_params)

cmap = LinearSegmentedColormap.from_list('interaction',
                                        ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])


def apa_analysis(apa):
    
    # remove outliers
    mean_arr = np.r_[[np.mean(arr) for arr in apa]]
    p99 = np.percentile(mean_arr, 99)
    p1 = np.percentile(mean_arr, 1)
    mask = (mean_arr < p99) & (mean_arr > p1)
    avg = apa[mask].mean(axis=0)
    lowerpart = avg[-5:,:5]
    ## z-score
    z = (avg[10, 10] - lowerpart.mean()) / lowerpart.std()
    
    return avg, z

queue = glob.glob('intermediate/*Pol2-control-loops.apa.*pkl')
vmin = 1
vmax = 1.5
for q in queue:
    outfig = os.path.split(q)[1].replace('.pkl', '.svg')
    if os.path.exists(outfig):
        continue

    apa = joblib.load(q)
    apa = np.r_[apa]
    avg, z = apa_analysis(apa)

    size = (0.88, 0.8)
    fig = plt.figure(figsize=size)
    width = 0.7
    Left = 0.1
    HB = 0.1
    HH = width * size[0] / size[1]
    ax = fig.add_axes([Left, HB, width, HH])
    sc = ax.imshow(avg, cmap=cmap, aspect='auto', interpolation = 'none',
                   vmax=vmax, vmin=vmin)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x = (xmin + xmax) / 2
    y = ymax + 1
    #ax.text(x, y, 'z-score: {0:.2g}'.format(z), va='top', ha='center', fontsize=5)
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(0.7)
    
    if 'Hi-C_FA_DSG_DpnII' in q:
        ax.set_title('Hi-C', fontsize=6, pad=3)
    elif 'MicroC.' in q:
        ax.set_title('Micro-C', fontsize=6, pad=3)
    elif 'ChIA-PET_CTCF.' in q:
        ax.set_title('CTCF ChIA-PET', fontsize=6, pad=3)
    elif 'ChIA-PET_Pol2.' in q:
        ax.set_title('Pol2 ChIA-PET', fontsize=6, pad=3)
    elif 'SPRITE.' in q:
        ax.set_title('DNA SPRITE', fontsize=6, pad=3)
    elif 'SPRITE_101_1000' in q:
        ax.set_title('DNA SPRITE (101-1000 fragments)', fontsize=6, pad=3)
    elif 'SPRITE_11_100' in q:
        ax.set_title('DNA SPRITE (11-100 fragments)', fontsize=6, pad=3)
    elif 'SPRITE_2_10' in q:
        ax.set_title('DNA SPRITE (2-10 fragments)', fontsize=6, pad=3)
    elif 'PLAC-Seq' in q:
        ax.set_title('H3K4me3 PLAC-Seq', fontsize=6, pad=3)
    else:
        ax.set_title('GAM', fontsize=6, pad=3)
    
    ax = fig.add_axes([Left+width+0.04, 0.72, 0.03, 0.15])
    fig.colorbar(sc, cax=ax, ticks=[vmin, vmax], format='%.3g')
    ax.tick_params(labelsize=5)

    plt.savefig(outfig, bbox_inches='tight', dpi=1000)
    plt.close()



