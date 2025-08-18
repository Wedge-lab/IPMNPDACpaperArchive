import pandas as pd, numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams

import shutup
shutup.please()

rcParams['font.family'] = 'Arial'

datapath = '/IPMNPDAC_WGS/Data/'

#updated15-4-23 sample label
precancerList =['case13_4','case10_2','case10_5','case12_S10','case12_S11',
                     'case12_S13','case12_S9','case15_1','case15_3','case2_S10',
                     'case3_1','case3_2','case3_4','case4_S1','case4_S3','case4_S4',
                     'case4_S5','case15_4','case16_2','case2_S2','case2_S4','case7_1']
cancerList = ['case6_S7','case7_4','case7_5','case9_S4','case4_S2','case11_S7',
                  'case11_S8','case13_3','case13_5','case15_10','case15_11','case16_4',
                   'case16_5','case3_5','case6_S8','case6_S9','case9_S2','case9_S3','case9_S6']

# snv data
snvdf =pd.read_csv(datapath + 'all41SNVTypeSampleCounts.csv')
snvprecancerDf = snvdf.query('samples==@precancerList')
snvprecancerDf.insert(1,'tumorStage', 'IPMN')
snvcancerDf = snvdf.query('samples==@cancerList')
snvcancerDf.insert(1, 'tumorStage','PDAC')
snvtumorStageDf = pd.concat([snvprecancerDf, snvcancerDf])
snvtumorStageDfb = snvtumorStageDf [list(snvtumorStageDf)[2:]]
snvtumorStageDfb.insert(0, 'tumorStage', list(snvtumorStageDf.tumorStage))
snvtumorStageDfb = snvtumorStageDfb.rename({'nonDriverCoding':'snvnonDriverCoding',
                                            'nonDriverRegulation':'snvnonDriverRegulation',
                                            'nonDriverIntronic':'snvnonDriverIntronic',
                                            'nonDriverintergenic':'snvnonDriverintergenic',
                                            'DriverNonCoding':'snvDriverNonCoding'}, axis=1)

snvtumorStageDfb = snvtumorStageDfb.reset_index(drop=True)

# indel data
indeldf =pd.read_csv(datapath + '41IndelTypeSampleCounts.csv')
indelprecancerDf = indeldf.query('samples==@precancerList')
indelprecancerDf.insert(1,'tumorStage', 'IPMN')
indelcancerDf = indeldf.query('samples==@cancerList')
indelcancerDf.insert(1, 'tumorStage','PDAC')
indeltumorStageDf = pd.concat([indelprecancerDf, indelcancerDf])
indeltumorStageDfb = indeltumorStageDf [list(indeltumorStageDf)[2:]]
indeltumorStageDfb.insert(0, 'tumorStage', list(indeltumorStageDf.tumorStage))
indeltumorStageDfb = indeltumorStageDfb.rename({'nonDriverCoding':'indelnonDriverCoding',
                                                'nonDriverRegulation':'indelnonDriverRegulation',
                                                'nonDriverIntronic':'indelnonDriverIntronic',
                                                'nonDriverintergenic':'indelnonDriverintergenic',
                                                'DriverNonCoding':'indelDriverNonCoding'}, axis=1)

indeltumorStageDfb = indeltumorStageDfb.reset_index(drop=True)
allDf = pd.concat([snvtumorStageDfb,indeltumorStageDfb], axis=1)
alldf = allDf.T.drop_duplicates().T

# plot
boxpps = dict(linestyle='-', linewidth=0, color='r')
medianpps = dict(linestyle='-', linewidth=1, color='r')

xt = alldf.boxplot(by='tumorStage',  medianprops=medianpps, sharey=False,
                   boxprops=boxpps,rot=0, grid=False, showfliers=False,
                   layout=(2,5), fontsize=10, return_type='both',figsize=(14,6),
                   patch_artist = True, column=list(alldf)[1:])

textposion = [100,680,2900,3300,43,0,46,250,310,0] #8 and 5 for no printing
pv=['p = 0.012','p = 0.005','p = 0.0038','p = 0.0034','p = 0.0004', '', 'p = 0.0045','p = 0.0047', 'p = 0.014', '']
colors = ['lightgreen',  'pink' ]
for i, (row_key, (ax, row)) in enumerate(xt.items()):
    ax.set_xlabel("")
    ax.set_title(row_key)
    ax.set_xticklabels("")
    ax.set_ylabel('Number of SVs')
    if i ==5 or i == 9:
       ax.text(0.65,textposion[i], '', fontsize=10) 
    else:
        ax.text(0.65,textposion[i], pv[i], fontsize=10)

    if row_key == 'snvnonDriverCoding':
        ax.set_ylabel('Number of SNVs')
    elif row_key == 'indelnonDriverCoding':
        ax.set_ylabel('Number of Indels')
    else:
        ax.set_ylabel('')
    for j,box in enumerate(row['boxes']):
        box.set_facecolor(colors[j])

plt.suptitle("")
plt.xticks([])

ipmn_patch = mpatches.Patch(color='lightgreen', label='IPMN')
pdac_patch = mpatches.Patch(color='pink', label='PDAC')
#plt.legend(handles=[ipmn_patch, pdac_patch], loc='best', frameon=False)
plt.tight_layout(pad=1)
plt.show()
