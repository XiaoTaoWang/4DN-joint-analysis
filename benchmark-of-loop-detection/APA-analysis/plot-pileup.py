import matplotlib, pickle, sys, glob, joblib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#ff9292','#ff6767','#F70000'])

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

def pileup(arrs):

    mean_arr = np.r_[[np.mean(arr) for arr in arrs]]
    p99 = np.percentile(mean_arr, 99)
    p1 = np.percentile(mean_arr, 1)
    mask = (mean_arr < p99) & (mean_arr > p1)
    pool = np.r_[arrs]
    avg = pool[mask].mean(axis=0)

    return avg
    
queue = glob.glob('*.loop-pileup.25kb.pkl')
vmin = 1
vmax = 1.8
for in_fil in queue:
    outf = in_fil.replace('.pkl', '.png')
    if not os.path.exists(os.path.join('selected', outf)):
        continue
    fig = plt.figure(figsize=(1.2,1.2))
    arrs = joblib.load(in_fil)
    avg = pileup(arrs)
    ax = fig.add_subplot(111)
    sc = ax.imshow(avg, vmin=vmin, vmax=vmax, interpolation='none', cmap=cmap)
    ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    c_ax = fig.add_axes([0.93, 0.7, 0.05, 0.18])
    cbar = fig.colorbar(sc, cax=c_ax, ticks=[vmin, vmax], format='%.3g')
    c_ax.tick_params(labelsize=5.5, pad=1, length=0.3)
    plt.savefig(outf, dpi=500, bbox_inches='tight')
    plt.close()

'''
for i, avg in enumerate([avg1, avg2, avg3]):
    outf = 'pileup-SVs.{0}.png'.format(i)
    norm = MidpointNormalize(vmin=avg.min(), vmax=avg.max(), midpoint=1)
    plt.imshow(avg, vmin=0.5, vmax=2, interpolation='none', cmap='seismic')
    plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.colorbar()
    plt.savefig(outf, dpi=300, bbox_inches='tight')
    plt.close()
'''
