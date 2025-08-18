import pandas as pd
from glob import glob

caseSampleClusterIDPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step4bcaseSampleClusterID/'
output = '/IPMNPDAC_WGS/Data/Data/sigDPC/step5bcaseClusterIDsum/'

dfs=[]
for fn in glob(caseSampleClusterIDPath + '*_chrPosCluster_ID.csv'):
    caseid = fn.split('/')[-1].split('_')[0]
    caseSampleClusterID_df = pd.read_csv(fn)
    caseSampleClusterID_df = caseSampleClusterID_df[['ID1',	'ID2',	'ID5',	'ID6',	'ID8',	'ID9',	'ID14', 'clusterNo']]
    sumIDCluster = caseSampleClusterID_df.groupby('clusterNo')[['ID1',	'ID2',	'ID5',	'ID6',	'ID8',	'ID9',	'ID14']].sum().round().reset_index()
    sumIDCluster.to_csv(output + '{}_clusterIDsum.csv'.format(caseid), index=0) 
    dfs.append(sumIDCluster)

