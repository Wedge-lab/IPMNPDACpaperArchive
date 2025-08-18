
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

svpath = '/IPMNPDAC_WGS/Data/'

svdf = pd.read_csv(svpath + '36treeSampleSV.csv')
multipleBranch = ['case10','case3','case16','case2','case9','case7','case6']
singleBranch = ['case4','case12','case13', 'case15']

mb_sv = svdf[svdf['samples'].str.contains('|'.join(multipleBranch))][["Total"]]
mb_sv.insert(0, 'samples', ['multipleBarnches']*mb_sv.shape[0])
sb_sv = svdf[svdf['samples'].str.contains('|'.join(singleBranch))][["Total"]]
sb_sv.insert(0, 'samples', ['singleBranch']*sb_sv.shape[0])
dfplot = pd.concat([mb_sv, sb_sv])

dfplot.to_csv(svpath + "SV_multiple_single_branches.csv", index=0)

boxpps = dict(linestyle='-', linewidth=0, color='r')
medianpps = dict(linestyle='-', linewidth=1, color='r')

ax = dfplot.boxplot(by='samples',  medianprops=medianpps, 
                    boxprops=boxpps,rot=0, grid=False, showfliers=False,
                    layout=(2,5), fontsize=10, return_type='both',figsize=(20,6),
                    patch_artist = True)

colors = ['lightgreen',  'pink' ]
for row_key, (ax, row) in ax.items():
    ax.set_title('')
    ax.set_ylabel("Total number of SVs")
    ax.set_xlabel("")
    ax.set_xticklabels(["Multiple branches", "Single branch"])
    #ax.set_yticklabels("")
    ax.text(1.7, 160, 'p = 0.0024', fontsize=10)
    for i,box in enumerate(row['boxes']):
        box.set_facecolor(colors[i])

plt.title('')
plt.suptitle('')
plt.xticks([] )
plt.yticks([])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ipmn_patch = mpatches.Patch(color='lightgreen', label='IPMN')
pdac_patch = mpatches.Patch(color='pink', label='PDAC')
plt.legend(handles=[ipmn_patch, pdac_patch], loc='best', frameon=False)

plt.show()
