import pandas as pd
from glob import glob

caseSampleClusterSBSPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step4caseSampleClusterSBS/'
output = '/IPMNPDAC_WGS/Data/Data/sigDPC/step5caseClusterSBSsum/'

dfs=[]
for fn in glob(caseSampleClusterSBSPath + '*_chrPosCluster_sbs.csv'):
    caseid = fn.split('/')[-1].split('_')[0]
    caseSampleClusterSBS_df = pd.read_csv(fn)
    caseSampleClusterSBS_df = caseSampleClusterSBS_df[['SBS1', 'SBS2',	'SBS5',	'SBS13','SBS17a', 
                                                       'SBS17b','SBS28', 'SBS40', 'clusterNo']]
    sumSBSCluster = caseSampleClusterSBS_df.groupby('clusterNo')[['SBS1', 'SBS2',	'SBS5',	'SBS13','SBS17a', 
                                                                  'SBS17b','SBS28', 'SBS40']].sum().round().reset_index()
    sumSBSCluster.to_csv(output + '{}_clusterSBSsum.csv'.format(caseid), index=0) 
    dfs.append(sumSBSCluster)

