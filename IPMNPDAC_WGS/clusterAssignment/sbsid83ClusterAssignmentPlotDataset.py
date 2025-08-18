import pandas as pd
from glob import glob

sbsClusterPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step5caseClusterSBSsum/'
idClusterPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step5bcaseClusterIDsum/'
output = '/IPMNPDAC_WGS/Data/Data/sigDPC/step6SBS_ID_cluster4plot/'

sbsCluster = glob(sbsClusterPath + '*_clusterSBSsum.csv')
idCluster = glob(idClusterPath + '*_clusterIDsum.csv')

for fsbs, fid in zip(sbsCluster, idCluster):
    caseid = fsbs.split('/')[-1].split('_')[0]
    sbscluster = pd.read_csv(fsbs)
    idcluster = pd.read_csv(fid)
    sbs_id_cluster = pd.merge(sbscluster, idcluster, on='clusterNo', how='left').fillna(0)
    sbs_id_cluster.to_csv(output + '{}_msDPC_SBS96_ID83.csv'.format(caseid), index=0)
