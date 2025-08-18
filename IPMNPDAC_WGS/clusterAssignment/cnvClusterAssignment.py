import os
import pandas as pd
import pyranges as pr

datapath = os.path.expanduser('/IPMNPDAC_WGS/Data/cnvDriverCluster/')

inputx = 'case16_pindelSNVmsDPC.csv'
outputx = 'case16_cnvDriverCluster.csv'

# CNV driver regions
[['hgnc_symbol','chromosome_name','start_position','end_position']]
CNVdriver = pd.read_csv(datapath + 'cnv_driversChrPos.csv')[['hgnc_symbol','chromosome_name','start_position','end_position']]
CNVdriver['chromosome_name'] = CNVdriver['chromosome_name'].astype(str)

# SNV clusters as intervals (start=end=pos)
clusterDf = pd.read_csv(datapath + inputx)[['chr', 'pos', 'most.likely.cluster']]
clusterDf['chr'] = clusterDf['chr'].astype(str)
clusterDf['start'] = clusterDf['pos']
clusterDf['end'] = clusterDf['pos'] + 1  # pyranges expects intervals

# Convert to PyRanges
gr_cnv = pr.PyRanges(CNVdriver.rename(columns={'chromosome_name':'Chromosome','start_position':'Start','end_position':'End'}))
gr_snv = pr.PyRanges(clusterDf.rename(columns={'chr':'Chromosome','start':'Start','end':'End'}))

# Find overlaps
overlaps = gr_snv.join(gr_cnv).df

# Add chr_pos field
overlaps['chr_pos'] = overlaps['Chromosome'].astype(str) + "_" + overlaps['pos'].astype(str)

# Save
overlaps.to_csv(datapath + outputx, index=False)
