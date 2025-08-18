import os
import glob
import shutil
import pandas as pd

datapath = '/mnt/d/1a_18_04_MU/5B_PancreaticCancer/4_IPMN_WGS/6b_callersDPC/10_CNVpindelSNVmsDPC/1_msDPC_mostUextract/'
outputpath = '/mnt/b/1_3July23VersionDataPlotCode/gitdata/Data/sigDPC/step3caseChrPosCluster/'

dfs=[]
for fd in (glob.glob(datapath + '*pindelSNVmsDPC')):
    caseID = fd.split('/')[-1].split('_')[0]
    snvFile = '{}__DP_and_cluster_info_0.01.txt'.format(caseID)#the file here could be snv+pindel
    clustFile = '{}__union_filtered_SNVs.txt'.format(caseID)
    pathx = os.path.join(os.getcwd(), fd)   
    snvFilex = os.path.join(pathx, snvFile)
    clustFilex = os.path.join(pathx, clustFile)
    dfchrpos = pd.read_csv(snvFilex, sep='\t')
    dfcluste = pd.read_csv(clustFilex, sep='\t')
    dfIDcluster = pd.concat([dfchrpos, dfcluste], axis=1)
    dfIDcluster.insert(0, 'chr_pos', dfIDcluster['chr'].astype(str) + '_' + dfIDcluster['pos'].astype(str))
    dfIDcluster =  dfIDcluster[['chr_pos', 'chr', 'pos', 'most.likely.cluster']]
    dfIDcluster = dfIDcluster.rename(columns = {'most.likely.cluster':'clusterNo'})
    dfs.append(dfIDcluster)
    dfIDcluster.to_csv(outputpath + '{}_chrPosCluster.csv'.format(caseID), index=0)

