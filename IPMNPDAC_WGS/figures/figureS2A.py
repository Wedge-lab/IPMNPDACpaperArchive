import pandas as pd
import matplotlib.pyplot as plt
from natsort import natsort_keygen
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

sfilepath = '/IPMNPDAC_WGS/Data/'

## 1a) input SBS96
sbs96all = pd.read_csv(sfilepath + 's41SBSsignatureCount.csv')
sbsplotdf = sbs96all[['typeTumorSample', 'SBS1','SBS2','SBS5',
                      'SBS13', 'SBS17a','SBS17b', 'SBS28', 'SBS40']]
sbsplotdf = sbsplotdf.sort_values(by="typeTumorSample", key=natsort_keygen())
sbsplotdf['typeTumorSample'] = [x[3:] for x in sbsplotdf.typeTumorSample]
sbsplotdfb = sbsplotdf.drop(['typeTumorSample'],axis=1)
sbsplotdf_total = sbsplotdfb.sum(axis=1)
sbsplotdf2rate = sbsplotdf[sbsplotdf.columns[1:]].div(sbsplotdf_total, 0) * 100
sbsplotdf2rate['typeTumorSample'] = sbsplotdf.typeTumorSample

# 1b) input ID83
id83all = pd.read_csv(sfilepath + 's41indelSignatureCount.csv')          
idplotdf = id83all[['typeTumorSample', 'ID1','ID2','ID5','ID6','ID8','ID9','ID14']]
idplotdf = idplotdf.sort_values(by="typeTumorSample", key=natsort_keygen())
idplotdf['typeTumorSample'] = [x[3:] for x in idplotdf.typeTumorSample]

idplotdfb = idplotdf.drop(['typeTumorSample'],axis=1)
idplotdf_total = idplotdfb.sum(axis=1)
idplotdf2rate = idplotdf[idplotdf.columns[1:]].div(idplotdf_total, 0) * 100
idplotdf2rate['typeTumorSample'] = idplotdf.typeTumorSample

# 1c input SV32
svall = pd.read_csv(sfilepath +'s41SVsignatureCount.csv')
svplotdf = svall[['typeTumorSample', 'SV2','SV4','SV5','SV7', 'SV9', 'SV10']]
svplotdf = svplotdf.sort_values(by="typeTumorSample", key=natsort_keygen())
svplotdf['typeTumorSample'] = [x[3:] for x in svplotdf.typeTumorSample]

svplotdfb = svplotdf.drop(['typeTumorSample'],axis=1)
svplotdf_total = svplotdfb.sum(axis=1)
svplotdf2rate = svplotdf[svplotdf.columns[1:]].div(svplotdf_total, 0) * 100
svplotdf2rate['typeTumorSample'] = svplotdf.typeTumorSample

# 1d) input CN48 
cnall = pd.read_csv(sfilepath +'s41CNsignature.csv')
cnplotdf = cnall[['typeTumorSample', 'CN1','CN9','CN24','CNV48B']]
cnplotdf = cnplotdf.sort_values(by="typeTumorSample", key=natsort_keygen())
cnplotdf['typeTumorSample'] = [x[3:] for x in cnplotdf.typeTumorSample]

cnplotdfb = cnplotdf.drop(['typeTumorSample'],axis=1)
cnplotdf_total = cnplotdfb.sum(axis=1)
cnplotdf2rate = cnplotdf[cnplotdf.columns[1:]].div(cnplotdf_total, 0) * 100
cnplotdf2rate['typeTumorSample'] = cnplotdf.typeTumorSample

# 2) set up plot
fig, axes = plt.subplots(1, 4, figsize=(20, 10), sharey=True)
colorSBS = ['blue', 'orange',  'green',  'red', 'purple', 'brown', 'pink', 'gray']
colorID = ['blue', 'orange', 'olive', 'cyan', 'lime', 'brown','fuchsia']
colorSV = ['blue', 'orange', 'olive', 'cyan', 'lime', 'brown']
colorCN = ['blue', 'green', 'gray', 'orange']

# 3a) plot sbs sigs
sbsplotdf2rate.plot(ax=axes[0], x = 'typeTumorSample', kind = 'barh', stacked = True,            
                    mark_right = True, legend=True, color=colorSBS, width=1.0)
axes[0].legend(bbox_to_anchor=(0.26, -0.1),fontsize=12)
axes[0].set_xlabel('Proportion of SBS96 in Each Sample', fontsize=14,weight='bold')

# 3b) plot Id sigs
idplotdf2rate.plot(ax=axes[1], x = 'typeTumorSample', kind = 'barh', stacked = True,            
                   mark_right = True, legend=True, color=colorID, width=1.0)

axes[1].legend(bbox_to_anchor=(0.22, -0.1),fontsize=12)
axes[1].set_xlabel('Proportion of ID83 in Each Sample',fontsize=14,weight='bold')

# 3c plot SV sigs
svplotdf2rate.plot(ax=axes[2], x = 'typeTumorSample', kind = 'barh', stacked = True,            
                   mark_right = True, legend=True, color=colorSV, width=1.0)
axes[2].set_xlabel('Proportion of SV32 in Each Sample',fontsize=14,weight='bold');
axes[2].legend(bbox_to_anchor=(0.23, -0.1),fontsize=12)

# 3d) plot CN sigs
cnplotdf2rate.plot(ax=axes[3], x = 'typeTumorSample', kind = 'barh', stacked = True,            
                 mark_right = True, legend=True, color=colorCN, width=1.0)
axes[3].set_xlabel('Proportion of CN48 in Each Sample',fontsize=14,weight='bold');
axes[3].legend(bbox_to_anchor=(0.26, -0.1), fontsize=12)

# 4) spacing
fig.tight_layout(pad=2.0)

plt.show();
