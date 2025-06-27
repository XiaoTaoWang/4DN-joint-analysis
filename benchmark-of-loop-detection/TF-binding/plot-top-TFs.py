import joblib, matplotlib
import numpy as np
import matplotlib.pyplot as plt

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

ctcf = 'without-CTCF'
data = joblib.load('{0}.TF-enrich.pkl'.format(ctcf))
colors = {
    4: '#33A02C',
    5: '#1F78B4',
    3: '#FDBF6F',
    2: '#B2DF8A'
}
for ci in [2, 3, 4, 5]:
    outfil = '{0}.C{1}.svg'.format(ctcf, ci+1)
    sort_table = []
    for tf in data[ci]:
        if tf in ['POLR2AphosphoS5']:
            continue
        t_per, enrich = data[ci][tf]
        sort_table.append([enrich, t_per, tf])
    
    sort_table.sort(reverse=True)

    # Select top 10 TFs
    top10 = sort_table[:12]
    enrich_scores = [x[0] for x in top10]
    t_per_scores = [x[1] for x in top10]
    tf_names = [x[2] for x in top10]

    fig = plt.figure(figsize=(1.7, 1))
    ax = fig.add_subplot(111)
    ax.bar(list(range(1, 13)), enrich_scores, color=colors[ci], ec='none', alpha=0.9)
    for x, y, text in zip(list(range(1, 13)), enrich_scores, t_per_scores):
        ax.text(x, y, '{0:.1%}'.format(text), va='bottom', ha='center', fontsize=5, rotation=270)
    
    ax.set_ylabel('Fold Enrichment', fontsize=6)
    ax.set_xlabel('TFs', fontsize=6)
    ax.set_title('Cluster {0}'.format(ci+1), fontsize=6)
    ax.set_xticks(list(range(1, 13)))
    ax.set_xticklabels(tf_names, rotation=45, ha='right', fontsize=5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)
    ax.xaxis.set_tick_params(width=0.8, pad=1.5, labelsize=5)
    ax.yaxis.set_tick_params(width=0.8, pad=1.5, labelsize=5)

    plt.savefig(outfil, dpi = 300, bbox_inches='tight')
    plt.close()