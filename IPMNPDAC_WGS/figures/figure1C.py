import os
import pandas as pd
import comuta 
#changed line 996 from new_cats = side_cats - paired_cats to new_cats = list(side_cats - paired_cats) 
#changed line 1002 from  missing_categories = paired_cats - side_cats to  missing_categories = list(paired_cats - side_cats)
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'  # Affects all text elements

driverheatmappath = '/IPMNPDAC_WGS/Data/'

# 1) Data processing
# 1a) read in snv-indel-brass-cnv data
drivermutsDf = pd.read_csv(driverheatmappath + '41sampleDriverAllmutsDirectHeatmap.csv')
drivermutsDf = drivermutsDf.query('value != "Silent"') 

# 1b) read in TMB data
tmbDf = pd.read_csv(driverheatmappath + 'tmbdf4driverheatmapPlot.csv')
sampleTumourTMBorderList = list(tmbDf ['sample'])

# 1.c) TMB data for top bar
tmbdf = pd.read_csv(driverheatmappath + '41sampleTMBdirectHeadmap.csv')

# 1d) read in purity data
cellularityDf = pd.read_csv(driverheatmappath + '41samplesCellularity4alldriverMutsHeatmap.csv')

# 1e) read in bar data
sidebarDf = pd.read_csv(driverheatmappath + '41sample_side_bar_data4alldrivermutsheatmap.csv')

category_order = list(sidebarDf.category)
sampleOrderForlabel = sampleTumourTMBorderList

###################################################################
mut_mapping = {'deletion_3UTRexon': 'lightgreen',
               'deletion_exon': 'deepskyblue',
               'deletion_nonExon': 'lime',
               'inversion_3UTRexon': 'darkkhaki',
               'inversion_nonExon': 'wheat',
               'tandem-duplication_nonExon':'thistle',
               'translocation_5UTRexon': 'turquoise',
               'translocation_nonExon': 'dodgerblue','Loss':'violet',
               'Invasive': 'maroon', 'Non-invasive': 'orange',
               'Pancreatitis': 'gray', 'IPMN_LGD': 'mediumslateblue',
               'IPMN_HGD': 'darkslateblue',
               'IPMN_HGD_PDAC': 'lightcoral', 'PDAC': 'maroon',
               'downstream': 'green', '3UTR': 'blue', 'nonsense': 'red',
               'ess_splice_Del': "purple", 'frameshift_Del': "teal",
               'frameshift_Ins': "slategray", 'upstream': 'lightblue',
               'missense': 'cyan',
               #'Silent': 'orange',
               '0.0_0.0': {'facecolor': 'grey', 'alpha': 0.05},
               'Absent': {'facecolor': 'grey', 'alpha': 0.05}}

##################################################################
# 2) plot
# 2a) plot heatmap
toy_comut = comuta.CoMut()
# add indicator, can add this after toy_comut.add_categorical_data, but plot on top
indicator_kwargs = {'color': 'red', 'marker': 'o',
                    'linewidth': 1, 'markersize': 5}

toy_comut.add_categorical_data(drivermutsDf, name='mut type',
                               mapping=mut_mapping,
                               category_order=category_order)

# 2b) plot sidebar
side_mapping = {'mutnb': 'pink'}
bar_kwargs = {'height': 0.8}
toy_comut.add_side_bar_data(sidebarDf, 
                            paired_name ='mut type', 
                            name='MutFq',
                            mapping =side_mapping, 
                            position ='right', bar_kwargs=bar_kwargs)

## 2c) plot cellularity
cat_mapping = {'Absent': {'facecolor': 'red'}}
value_range = (0, 1)
toy_comut.add_continuous_data(cellularityDf, 
                              name = 'Purity', 
                              mapping = 'gray_r', 
                              cat_mapping = cat_mapping, value_range = value_range)

# 2d) plot tmb bar
bar_mapping = {'Nonsynonymous': 'green', 'Synonymous': 'pink'}
bar_kwargs = {'width': 0.8, 'edgecolor': 'green'}

toy_comut.add_bar_data(tmbDf, name = 'Mutation burden', mapping = bar_mapping,  bar_kwargs = bar_kwargs)
                       #ylabel = 'Muts/Mb')

toy_comut.plot_comut(figsize=(16, 12), x_padding=0.04,hspace =0.03, 
                     wspace = 0.03, y_padding=0.04, tri_padding=0.03)

toy_comut.axes['Mutation burden'].set_ylabel('TMB', rotation = 'horizontal', 
                                             ha = 'right', va = 'center', y = 0.3)
toy_comut.axes['Mutation burden'].set_ylim(ymin=0, ymax=10)
toy_comut.axes['Mutation burden'].axvline(0, color='black')
#toy_comut.axes['mut type'].set_xlabel('xt', rotation = 45)
toy_comut.axes['mut type'].set_xticklabels(labels=sampleOrderForlabel,rotation = 45,horizontalalignment='right')

custom_rcParams = {'font.size': 12}
rcParams.update(custom_rcParams)
toy_comut.add_unified_legend(bbox_to_anchor = (1.8, -4))
#toy_comut.add_unified_legend()
plt.subplots_adjust(bottom=0.001, right=0.9, top=0.99,left=0.02);
plt.savefig(driverheatmappath + 'driverHeatMap.png', dpi=300, bbox_inches='tight')
