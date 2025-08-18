import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams

import shutup
shutup.please()

rcParams['font.family'] = 'Arial'

## TMB bar plot data
tmbDf = pd.read_csv('/IPMNPDAC_WGS/Data/tmbplotPermutest.csv')[['tumorStage',  'snv_indel_TMB']]
ipmn = tmbDf.query('tumorStage=="IPMN"')
pdac = tmbDf.query('tumorStage=="PDAC"')

## boxplot
data = [ipmn.snv_indel_TMB, pdac.snv_indel_TMB]
fig = plt.figure(figsize =(5,5))
ax = fig.add_subplot(111)

bp = ax.boxplot(data, patch_artist = True,
                notch ='True', vert =90)
colors = ['lightgreen', 'pink']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='blue',
                linewidth = 1.5,
                linestyle =":")
# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='purple',
            linewidth = 2)
# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='red',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='black',
              alpha = 0.5)   

#plt.title("TMB profile of IPMN-PDAC", fontsize=14, weight='bold')
plt.ylabel('Tumor mutation burden', fontsize=14)
# Removing top axes and right axes
# ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# show plot
import matplotlib.patches as mpatches
ipmn_patch = mpatches.Patch(color='lightgreen', label='IPMN')
pdac_patch = mpatches.Patch(color='pink', label='PDAC')
plt.legend(handles=[ipmn_patch, pdac_patch], loc='best', frameon=False)
plt.xticks([])
plt.show()
