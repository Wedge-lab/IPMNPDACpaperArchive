import os
import numpy as np
import pandas as pd
import comuta
from natsort import index_natsorted,natsorted,natsort_keygen
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

import shutup
shutup.please()

datapath = '/IPMNPDAC_WGS/Data/'

# Data processing
# 1.a) sampleID
sampleDf = pd.read_csv(datapath+ 'sample_tumorType.csv')
sample_orderTumor = list(sampleDf.typeTumosample)
sampleGroup = list(sampleDf.group)

sampleOldName = list(sampleDf.samples)
sampleNewName = list(sampleDf.typeTumosample)

sampleIndicator = pd.DataFrame({'sample':sample_orderTumor, 'group':sampleGroup})

sampleIndicatorSampleID = sampleIndicator
sampleIndicatorSampleID['sampleID'] = [a.split('_')[-2]+'_'+a.split('_')[-1] for a in list(sampleIndicatorSampleID['sample'])]
sampleIndicatorSampleID = sampleIndicatorSampleID.rename(columns={'sample':'sampleTumor', 'sampleID':'sample'})

# 1.b) TMB data for top bar
tmbdf = pd.read_csv(datapath + '41sampleTMBdirectHeadmap.csv')
sampleTumorTMB = sampleIndicatorSampleID.merge(tmbdf)

tumourGrade=['t1','t2','t3','t4','t5']
sampleTumorTMBOrder=[]
for td in tumourGrade:
    dfx = sampleTumorTMB[sampleTumorTMB['sampleTumor'].str.contains(td)]
    dfx = dfx.sort_values(by='Nonsynonymous',ascending=True)
    sampleTumorTMBOrder.append(dfx)

sampleTumorTMBOrderDf = pd.concat(sampleTumorTMBOrder)
sampleTumourTMBorderList = list(sampleTumorTMBOrderDf['sample'])
tmbdf2 = sampleTumorTMBOrderDf[['sample', 'Nonsynonymous',  'Synonymous']]
tmbdf2.to_csv(datapath + 'tmbdf4driverheatmapPlot.csv', index=0)

# 1.c) CNV data
totalDf = pd.read_csv(datapath + 'ipmn_gisticData.txt', sep='\t')
cytobandOnly = [a.split(':')[1] for a in totalDf.Cytoband]
heatmapDf = totalDf[['Tumor_Sample_Barcode','Variant_Classification']]
heatmapDf.insert(2, 'Cytoband', cytobandOnly)
heatmapDf = heatmapDf.rename({'Tumor_Sample_Barcode':'sample', 
                              'Cytoband':'category', 
                              'Variant_Classification':'value'}, axis=1) 

heatmapDf = heatmapDf.replace(sampleOldName,sampleNewName)
heatmapDf = heatmapDf.drop_duplicates()

# 1.d) read in sample type data
tumourType = pd.read_csv(datapath + 'all41TumourType.csv')

# 1.e) fq side bar  
heatmapDf = pd.concat([heatmapDf,tumourType])
heatmapDf = heatmapDf.sort_values(by=["sample"], key=natsort_keygen())

# 1.f) cellularity top grad bar
purity_data = pd.read_csv(datapath + '41samplesCellularity4alldriverMutsHeatmap.csv')
#purity_data.insert(1, 'category', 'Purity')

# 2) plot
# 2.a) sample order and label
sampleOrderForlabel = sampleTumourTMBorderList 
heatmapDf['sample'] = [y.split('_')[-2]+'_'+y.split('_')[-1] for y in heatmapDf['sample']]

heatmapDf['sample'] = pd.Categorical(heatmapDf ['sample'], ordered=True, categories=sampleTumourTMBorderList)
heatmapDf = heatmapDf.sort_values('sample')
heatmapDf2 = heatmapDf.reset_index()
heatmapDf2 = heatmapDf2[['sample', 'value', 'category']]
side_bar_data = pd.DataFrame({'count':heatmapDf2.groupby( ["category"] ).size()}).reset_index()
side_bar_data['-log(Q)'] = side_bar_data['count'].div(0.41).round(1)
side_bar_data = side_bar_data.drop('count', axis=1)
side_bar_data.loc[side_bar_data['category'] == 'Pathology', ['count',  '-log(Q)']] = 0
side_bar_data.loc[side_bar_data['category'] == 'Invasive', ['count',  '-log(Q)']] = 0
side_bar_data = side_bar_data.sort_values(by=['-log(Q)'], 
                                          key=lambda x: np.argsort(index_natsorted(side_bar_data['-log(Q)'])))
category_order = side_bar_data.category
value_order = ['Del','Amp']

purity_data = purity_data.set_index('sample')
purity_data = purity_data.reindex(index = sampleOrderForlabel)
purity_data = purity_data.reset_index()
purity_data  = purity_data.rename(columns = {'cellularity':'value'}) #must keep the same column name reqored by the package

# 2.b) plot colour
mut_mapping = {'Amp':'skyblue', 'Del':'violet', 
               'Absent':{'facecolor':'grey', 'alpha':0.05},
               'Invasive':'maroon', 'Non-invasive':'orange',
               'Pancreatitis':'gray','IPMN_LGD':'mediumslateblue', 
               'IPMN_HGD':'darkslateblue', 
               'IPMN_HGD_PDAC':'lightcoral', 'PDAC':'maroon'}

side_mapping = {'-log(Q)':'pink'}
sidebar_kwargs = {'height': 0.8}

bar_mapping = {'Nonsynonymous': 'green', 'Synonymous': 'pink'}
bar_kwargs = {'width': 0.8, 'edgecolor': 'green'}

cat_mapping = {'Absent': {'facecolor': 'red'}}
value_range = (0, 1)

# 2.c) add heat main map
toy_comut = comuta.CoMut()
toy_comut.add_categorical_data(heatmapDf2, name='Mutation type', 
                               mapping=mut_mapping,
                               category_order=category_order,
                               value_order = value_order)

# 2.d) add side bar
toy_comut.add_side_bar_data(side_bar_data, paired_name ='Mutation type', name ='MutFq',
                            mapping = side_mapping, xlabel = 'Fq (%)', 
                            position = 'right', bar_kwargs = sidebar_kwargs)
# 2.e) plot cellularity
cat_mapping = {'Absent': {'facecolor': 'red'}}
value_range = (0, 1)
toy_comut.add_continuous_data(purity_data, 
                              name = 'Purity', 
                              mapping = 'gray_r', 
                              cat_mapping = cat_mapping, value_range = value_range)
# 2.f) add top bar
toy_comut.add_bar_data(tmbdf2, name = 'Mutation burden', 
                       mapping = bar_mapping,  
                       bar_kwargs = bar_kwargs,
                       ylabel ='TMB')

# 2.g) general adjust
toy_comut.plot_comut(figsize=(12, 10), x_padding=0.04,hspace =0.03, 
                     wspace = 0.03, y_padding=0.04, tri_padding=0.03)

toy_comut.axes['Mutation type'].set_xticklabels(labels=sampleOrderForlabel,rotation = 45,horizontalalignment='right')
toy_comut.axes['Mutation burden'].set_ylabel('TMB', rotation = 'horizontal', 
                                             ha = 'right', va = 'center', y = 0.3)
toy_comut.axes['Mutation burden'].set_ylim(ymin=0, ymax=10)
toy_comut.axes['Mutation burden'].axvline(0, color='black')

custom_rcParams = {'font.size': 10}
rcParams.update(custom_rcParams)

toy_comut.add_unified_legend(bbox_to_anchor = (1.5, -4))

plt.subplots_adjust(bottom=0.001, right=0.9, top=0.999,left=0.05);
