import os, glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker as mticker
from scipy.stats import mannwhitneyu

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

def get_loop_sizes(fils):

    sizes = []
    for fil in fils:
        with open(fil, 'r') as source:
            for line in source:
                c1, s1, e1, c2, s2, e2 = line.rstrip().split()[:6]
                p1 = (int(s1) + int(e1)) // 2
                p2 = (int(s2) + int(e2)) // 2
                sizes.append(abs(p2 - p1))
    
    return sizes

EP_loop_sizes = get_loop_sizes(['H1ESC.EP-loops.bedpe', 'HFFc6.EP-loops.bedpe'])
CTCF_loop_sizes = get_loop_sizes(['H1ESC.CTCF-loops.bedpe', 'HFFc6.CTCF-loops.bedpe'])
print('median size of EP loops:', np.median(EP_loop_sizes))
print('median size of CTCF loops:', np.median(CTCF_loop_sizes))
print(mannwhitneyu(EP_loop_sizes, CTCF_loop_sizes))

fig = plt.figure(figsize=(1, 1.2))
ax = fig.add_subplot(111)
EP_loop_sizes = np.log10(np.r_[EP_loop_sizes])
CTCF_loop_sizes = np.log10(np.r_[CTCF_loop_sizes])
colors = ['#FF7F00', '#33A02C', '#1F78B4']

bp = sns.violinplot(data=[EP_loop_sizes, CTCF_loop_sizes], ax=ax, palette=colors[1:])
for violin in ax.collections[::2]:
    violin.set_alpha(0.7)
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
#ax.yaxis.set_ticks([np.log10(x) for p in range(1,8) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Genomic distance between loop anchors (bp)', fontsize=6)
#ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['EP loops', 'CTCF loops'], fontsize=6, rotation=30, ha='right')
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.xaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.yaxis.set_tick_params(width=1, labelsize=6, pad=2)
ax.set_axisbelow(True)
plt.savefig('size-distributions-violin.svg', dpi=500, bbox_inches='tight')
plt.close()