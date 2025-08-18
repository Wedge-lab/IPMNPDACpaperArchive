import pandas as pd
from glob import glob

seqinfoSBSprobPath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step2seqInfoSigOut/'
caseChrPosClusterpath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step3caseChrPosCluster/'
outpath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step4caseSampleClusterSBS/'

# 1) process seqinfoSBSprob by adding column case_chrs_pos
seqinfoSBSprob = pd.read_csv(seqinfoSBSprobPath + '41seqInfoSBSprob.csv',low_memory=False)
seqinfoSBSprob = seqinfoSBSprob[['samples','chrs','pos', 'SBS1','SBS2','SBS5','SBS13','SBS17a','SBS17b','SBS28','SBS40']]
seqinfoSBSprob.insert(1, 'caseid', [x.split('_')[0] for x in seqinfoSBSprob.samples])
seqinfoSBSprob.insert(0, 'case_chrs_pos', seqinfoSBSprob['caseid'].astype(str)+'_'+ seqinfoSBSprob['chrs'].astype(str)+'_'+seqinfoSBSprob['pos'].astype(str))
seqinfoSBSprob = seqinfoSBSprob[['case_chrs_pos','samples','SBS1','SBS2','SBS5','SBS13','SBS17a','SBS17b','SBS28','SBS40']]

# 2) process caseChrPosCluster by adding column case_chrs_pos, and merge

dfs = []
for fn in glob(caseChrPosClusterpath + '*_chrPosCluster.csv'):
    caseid = fn.split('/')[-1].split('_')[0]
    chrPosCluster_df = pd.read_csv(fn)
    chrPosCluster_df.insert(0, 'case_chrs_pos',  caseid + '_' + chrPosCluster_df['chr_pos'].astype(str))
    chrPosCluster_df = chrPosCluster_df[['case_chrs_pos', 'clusterNo']]
    chrPosCluster_sbs = seqinfoSBSprob.merge(chrPosCluster_df)
    chrPosCluster_sbs = chrPosCluster_sbs.drop_duplicates(subset='case_chrs_pos', keep='first')
    chrPosCluster_sbs.to_csv(outpath + '{}_chrPosCluster_sbs.csv'.format(caseid), index=0)
    dfs.append(chrPosCluster_sbs)
