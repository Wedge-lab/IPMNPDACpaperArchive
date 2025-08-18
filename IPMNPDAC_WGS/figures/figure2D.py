import os
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt

#pd.options.mode.chained_assignment = None
from matplotlib import rcParams
plt.rcParams['figure.dpi'] = 300
rcParams['font.family'] = 'Arial'

os.chdir('/IPMNPDAC_WGS/Data/sigDPC/step6SBS_ID_cluster4plot/')

#Helper function to display absolute numbers
def absolute_number(pct, all_vals):
    total = sum(all_vals)
    absolute = int(round(pct * total / 100.0))
    return f"{absolute}" if absolute > 0 else ''


for fn in glob('*_msDPC_SBS96_ID83.csv'):
    sbs_id = pd.read_csv(fn)
    caseID = fn.split('_')[0]
    clusterx = list(sbs_id.clusterNo)
    sig_labels = list(sbs_id)[1:]
    labels_outer = sig_labels[:7]
    labels_inner = sig_labels[7:]


    fig, ax = plt.subplots(1,sbs_id.shape[0], figsize=(1.7*sbs_id.shape[0], 1.7*sbs_id.shape[0]))
    #plt.subplots_adjust(wspace = 0.00, hspace= 0.2, bottom=0.12, right=0.99, top=0.99,left=0.005)

    ds=[]
    for i in range(sbs_id.shape[0]):
        x=sbs_id.iloc[[i]] 
        id83 = sbs_id.iloc[[i]].filter(regex=r'(ID[0-9])')
        sbs96 = sbs_id.iloc[[i]].filter(regex=r'(SBS[0-9])')
        testid = np.array(id83.iloc[0].values)
        testsbs = np.array(sbs96.iloc[0].values)
        ds.append(testid)

        size = 0.3

        outer_colors = ['limegreen','cyan','teal', 'greenyellow',
                    'olive', 'blue', 'blueviolet','deepskyblue']
        inner_colors = ['red','#FF00FF', '#C20078','#DDA0DD',
                        'lightpink', 'orange', 'silver']

        if testsbs.sum() > 0 :

            ax[i].pie(testsbs, radius=1, colors=outer_colors,
                      wedgeprops=dict(width=size, edgecolor='w'),
                      autopct=lambda pct: absolute_number(pct, testsbs.tolist()),
                      pctdistance=1.15,
                      textprops=dict(color='black', fontsize=6)) # number size
            if testid.sum() > 0:
                ax[i].pie(testid, radius=1-size, colors=inner_colors,
                          wedgeprops=dict(width=size, edgecolor='w'),
                          autopct=lambda pct: absolute_number(pct, testid.tolist()),
                          pctdistance=0.75,
                          textprops=dict(color='black', fontsize=6))
        else:
            if testid.sum() > 0:
                ax[i].pie(testid, radius=1-size, colors=inner_colors,
                          wedgeprops=dict(width=size, edgecolor='w'),
                          autopct=lambda pct: absolute_number(pct, testid.tolist()),
                          pctdistance=0.75,
                          textprops=dict(color='black', fontsize=6))

        #ax[i].set(title='Input cluster{}'.format(clusterx[i]))    
        ax[i].set_title('Cluster{}'.format(clusterx[i]),fontsize=8)
    allColor = outer_colors + inner_colors
    name_to_color = {name:color for name, color in zip(sig_labels, allColor)}
    handles = [plt.Rectangle((0, 0), 0, 0, color=name_to_color[name], label=name) for name in name_to_color]

    plt.legend(handles=handles, loc='upper center',  bbox_to_anchor=(-3.2, -0.05),
               fancybox=True, shadow=True, ncol=8, frameon=False, fontsize=8)

    #ax.set(aspect="equal", title='Pie plot with `ax.pie`')
    plt.text(-18, 2.5, '{}: Proportion of SBS and ID signatures in subclones'.format(caseID), fontsize = 12)
    plt.show();

