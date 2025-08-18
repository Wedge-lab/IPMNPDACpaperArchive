import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

treeSample = ['case6_S7','case6_S9',
              'case7_S1','case7_S4','case7_S5',
              'case9_S2','case9_S3','case9_S4',
              'case10_S2','case10_S5',
              'case12_S1','case12_S2','case12_S3','case12_S5',
              'case13_S3','case13_S4','case13_S5',
              'case16_S2','case16_S4',
              'case2_S2','case2_S4','case2_S10',
              'case4_S1','case4_S2','case4_S3','case4_S4', 'case4_S5',
              'case15_S1','case15_S3', 'case15_S4','case15_S10','case15_S11',]

colorx =['lightgreen']*19 +['pink']*13
# 1) import dataset
os.chdir(r'/mnt/b/1_3July23VersionDataPlotCode/NC_revision/pathwayCelltype')
cellTypes = pd.read_csv('ESTIMATE_EPIC_scoresb.csv')
treeCells = cellTypes.query('Samples==@treeSample')
treeCells = treeCells.rename(columns = {'B_cell':'B cell', 'Cancer associated_fibroblast':'Fibroblast',
                                        'T_cell_CD4+' : 'CD4+','T_cell_CD8+':'CD8+','Endothelial_cell': 'Endothelial cell',
                                        'NK_cell':'NK cell','uncharacterized_cell':'Uncharacterized'})
treeCells = treeCells.set_index('Samples')
treeCells = treeCells.reindex(treeSample) # Apply custom sort

cellnames = list(treeCells)
lbs = list(treeCells.index)
binx = np.arange(len(lbs))

# 2) Create a figure and axis objects
fig, axes = plt.subplots(nrows=len(cellnames), ncols=1, figsize=(12, 20))

# 3) Iterate over the data and create subplots with bar charts
for i, ax in enumerate(axes):
    # Get the data for the current subplot
    subplot_df = treeCells[[cellnames[i]]]

    # Set the title for the subplot
    #ax.set_title(f"Subplot {i+1}")

    # Create the bar chart for the current subplot
    ax.bar(subplot_df.index, subplot_df[cellnames[i]]*100, width = 0.3, color=colorx)
    ax.set_xticklabels([])
    ax.set_xticks([]) 
    ax.set_xlim(-0.5, binx.size-0.5)
    # Set the labels for the x-axis and y-axis
    #ax.set_xlabel('X-axis')
    ax.set_ylabel(cellnames[i], fontsize=14)

ax.set_xticks(range(len(lbs)))    
ax.set_xticklabels(lbs, rotation=45, ha='right')    
#plt.xlim([-0.5,binx.size-0.5])
# Adjust the spacing between subplots

#plt.tight_layout()

# Show the plot
plt.show()
