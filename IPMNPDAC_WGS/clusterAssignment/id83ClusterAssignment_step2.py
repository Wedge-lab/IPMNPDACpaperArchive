import pandas as pd
from glob import glob

seqinfoIDprobPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step2seqInfoSigOut/'
caseChrPosClusterpath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step3caseChrPosCluster/'
outpath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step4bcaseSampleClusterID/'

# 1) process seqinfoSBSprob by adding column case_chrs_pos
seqinfoIDprob = pd.read_csv(seqinfoIDprobPath + '41seqInfoIDprob.csv')

# 2) process caseChrPosCluster by adding column case_chrs_pos, and merge

dfs = []
for fn in glob(caseChrPosClusterpath + '*_chrPosCluster.csv'):
    caseid = fn.split('/')[-1].split('_')[0]
    chrPosCluster_df = pd.read_csv(fn)
    chrPosCluster_df.insert(0, 'case_chrs_pos',  caseid + '_' + chrPosCluster_df['chr_pos'].astype(str))
    chrPosCluster_df = chrPosCluster_df[['case_chrs_pos', 'clusterNo']]
    chrPosCluster_ID = seqinfoIDprob.merge(chrPosCluster_df)
    chrPosCluster_ID = chrPosCluster_ID.drop_duplicates(subset='case_chrs_pos', keep='first')
    chrPosCluster_ID.to_csv(outpath + '{}_chrPosCluster_ID.csv'.format(caseid), index=0)
    dfs.append(chrPosCluster_ID)
