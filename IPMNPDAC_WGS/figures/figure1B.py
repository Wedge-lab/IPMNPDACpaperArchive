
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'  # Affects all text elements

import shutup
shutup.please()

datapath = '/IPMNPDAC_WGS/Data/'

## 0) preparation samplelegend
dfSample = pd.read_csv(datapath + 'tumorTypetemplate.csv')

## 1) preparation SNV
plotSNV = pd.read_csv(datapath + '41samplesSNVNumberBarPlot.csv')

## 2) preparation indels  
plotIndel = pd.read_csv(datapath + 'total41NumberIntersectIndels.csv')
plotIndel['typeTumorSample'] = [x[3:] for x in plotIndel.typeTumosample]

## 3) preparation SVs 
plotBrass = pd.read_csv(datapath + 'total41NumberBrass.csv')
plotBrass = plotBrass.rename(columns = {'tandem-duplication':'tand-dupl', 'translocation':'trans'})
plotBrass['typeTumorSample'] = [x[3:] for x in plotBrass.typeTumorSample]

## 4) SBS96 preparation  
plotsbsig = pd.read_csv(datapath + '41sbsigReorder.csv')
plotsbsig['typeTumorSample'] = [x[3:] for x in plotsbsig.typeTumorSample]

## 5) ID signature   
plotidsig = pd.read_csv(datapath + '41idsigReorder.csv')
plotidsig['typeTumorSample'] = [x[3:] for x in plotidsig.typeTumorSample]

## 6) SV signature
plotsvsig = pd.read_csv(datapath +'41samplesSVsigNbBarPlot.csv')

## 7) CN signature read in and preparation 
plotCNsig = pd.read_csv(datapath + '41cnsigReorder.csv')
plotCNsig['typeTumorSample'] = [x[3:] for x in plotCNsig.typeTumorSample]

## 8) plot 
fig, axes = plt.subplots(1, 8, figsize=(30, 15), dpi=144, 
                         gridspec_kw={'width_ratios': [0.2,1,1,1,1,1,1,1]},sharey=True) #adjuct te width
plt.subplots_adjust(wspace = 0.15,bottom=0.15, right=0.98, top=0.97,left=0.07) #hspace

dfSample.plot(ax=axes[0],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               color={'Invasive': 'maroon', 'Noninvasive': 'orange',
               'Pancreatitis': 'gray', 'IPMN_LGD': 'mediumslateblue',
               'IPMN_HGD': 'darkslateblue','Invasivex':'white', 'Pathology':'white',
               'IPMN_HGD_PDAC': 'lightcoral', 'PDAC': 'maroon'},
                width=1.0)

axes[0].set_ylabel('Case tumour type',fontsize=16,weight='bold')
axes[0].tick_params(axis='y', labelsize=14)
axes[0].set_xlim(0, 3)
axes[0].set_xticklabels(['Invasive', 'Pathology'], rotation=90, horizontalalignment='right',fontsize=16)
axes[0].set_frame_on(False)
axes[0].legend(loc="upper center")
handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, loc='center', bbox_to_anchor=(15,-0.12), ncol=9,frameon=False,fontsize=16)
#############
plotSNV.plot(ax=axes[1], 
             kind = 'barh',
             color ='skyblue',
             edgecolor = "black",
             width = 1.0)
axes[1].set_xlabel('a. Number of SNVs',fontsize=16,weight='bold')
axes[1].tick_params(axis='x', labelsize=14)
axes[1].set_xlim(0, 8000)
axes[1].legend(fontsize=14, frameon=False)

#########################
plotIndel.plot(ax=axes[2],
            x = 'typeTumosample',
            kind = 'barh', 
            stacked = True,            
            mark_right = True,
            edgecolor = "black",
            width=1.0)
axes[2].set_xlabel('b. Number of Ins/Dels',fontsize=16,weight='bold')
axes[2].tick_params(axis='x', labelsize=14)
axes[2].set_xlim(0, 800)
axes[2].legend(fontsize=14,frameon=False)
#########################
plotBrass.plot(ax=axes[3],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               edgecolor = "black",
               width=1.0)
axes[3].set_xlabel('c. Number of SVs',fontsize=16, weight='bold')
axes[3].tick_params(axis='x', labelsize=14)
axes[3].set_xlim(0, 250)
axes[3].legend(fontsize=14,frameon=False)
########################
plotsbsig.plot(ax=axes[4],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               edgecolor = "black",
               width=1.0)
axes[4].set_xlabel('d. Number of SBS96',fontsize=16, weight='bold')
axes[4].tick_params(axis='x', labelsize=14)
axes[4].set_xlim(0, 8000)
axes[4].legend(fontsize=14, frameon=False)
#######################
plotidsig.plot(ax=axes[5],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               edgecolor = "black",
               width=1.0)
axes[5].set_xlabel('e. Number of ID83',fontsize=16, weight='bold')
axes[5].tick_params(axis='x', labelsize=14)
axes[5].set_xlim(0, 800)
axes[5].legend(fontsize=14, frameon=False)
#########################
plotsvsig.plot(ax=axes[6],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               edgecolor = "black",
               width=1.0)
axes[6].set_xlabel('f. Number of SV32',fontsize=16, weight='bold')
axes[6].tick_params(axis='x', labelsize=14)
axes[6].set_xlim(0, 350)
axes[6].legend(fontsize=14, frameon=False)

#######################
plotCNsig.plot(ax=axes[7],
               x = 'typeTumorSample',
               kind = 'barh', 
               stacked = True,            
               mark_right = True,
               edgecolor = "black",
               width=1.0)
axes[7].set_xlabel('g. Number of CN48',fontsize=16, weight='bold')
axes[7].tick_params(axis='x', labelsize=14)
axes[7].set_xlim(0, 350)
axes[7].legend(fontsize=14, frameon=False)

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
plt.show();
