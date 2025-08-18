import os
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

from matplotlib import rcParams

import shutup
shutup.please()

rcParams['font.family'] = 'Arial'

svpath = '/IPMNPDAC_WGS/Data/'

#updated15-4-23 sample label
precancerList =['case13_4','case10_2','case10_5','case12_S10','case12_S11',
                     'case12_S13','case12_S9','case15_1','case15_3','case2_S10',
                     'case3_1','case3_2','case3_4','case4_S1','case4_S3','case4_S4',
                     'case4_S5','case15_4','case16_2','case2_S2','case2_S4','case7_1']
cancerList = ['case6_S7','case7_4','case7_5','case9_S4','case4_S2','case11_S7',
                  'case11_S8','case13_3','case13_5','case15_10','case15_11','case16_4',
                   'case16_5','case3_5','case6_S8','case6_S9','case9_S2','case9_S3','case9_S6']

# input
brassdf = pd.read_csv(svpath + '41BrassTypeSampleCounts.csv')

brassprecancerDf = brassdf.query('samples==@precancerList')
brassprecancerDf.insert(1,'tumorStage', 'IPMN')
brasscancerDf = brassdf.query('samples==@cancerList')
brasscancerDf.insert(1, 'tumorStage','PDAC')
brasstumorStageDf = pd.concat([brassprecancerDf, brasscancerDf])
brasstumorStageDfb = brasstumorStageDf[list(brasstumorStageDf)[1:]]

# plot
boxpps = dict(linestyle='-', linewidth=0, color='r')
medianpps = dict(linestyle='-', linewidth=1, color='r')

xt = brasstumorStageDfb.boxplot(by='tumorStage',  medianprops=medianpps, sharey=True,
                                boxprops=boxpps,rot=0, grid=False, showfliers=False,
                                layout=(1,4), fontsize=10, return_type='both', figsize=(16,4),
                                patch_artist = True, column=list(brasstumorStageDfb)[1:])

colors = ['lightgreen',  'pink' ]
pv=['p = 0.1051','p = 0.0002','p = 0.0004','p = 0.0004']
colors = ['lightgreen',  'pink' ]
for i, (row_key, (ax, row)) in enumerate(xt.items()):
    ax.set_xlabel("")
    #ax.set_title(row_key)
    ax.set_xticklabels("")
    ax.set_ylabel('Number of SVs')
    ax.text(0.65, 65, pv[i], fontsize=12)
    for i,box in enumerate(row['boxes']):
        box.set_facecolor(colors[i])

plt.suptitle("")
ipmn_patch = mpatches.Patch(color='lightgreen', label='IPMN')
pdac_patch = mpatches.Patch(color='pink', label='PDAC')
plt.subplots_adjust(wspace=0.12, hspace=0.5)
#plt.legend(handles=[ipmn_patch, pdac_patch], loc='best', frameon=False)

plt.show()
