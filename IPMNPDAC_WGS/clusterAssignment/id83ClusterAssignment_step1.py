import os
from glob import glob
import pandas as pd
import numpy as np
from natsort import natsorted
from natsort import index_natsorted

sampleID = pd.read_csv('/mnt/b/1_3July23VersionDataPlotCode/Figure1/sampleID2025.csv')
sampleNewID = list(sampleID.ID2025)
sampleOldID = list(sampleID.previousID)

# 1) ID83_seqinfo
id83seqinfoPath = '/mnt/e/5_signatureResult34_41samples/42ID83/output/vcf_files/ID/'

idDFs = []
for fm in glob(id83seqinfoPath +'*seqinfo.txt'):
    idDf = pd.read_csv(fm, sep='\t', header=None)
    idDf.columns = ['samples', 'chrs', 'pos', 'mutationType', 'REF', 'ALT', 'orientation']
    idDf = idDf[['samples', 'chrs', 'pos', 'mutationType', 'orientation']]
    idMutTypeSeq = [f[2:] for f in idDf.mutationType]
    idDf.insert(4, 'mutationTypeSeq', idMutTypeSeq)
    idDFs.append(idDf)   

IDdf = pd.concat(idDFs)  
IDdf_sampleMut = [str(g)+ '_' + str(h) for g, h in zip(IDdf['samples'], IDdf['mutationTypeSeq'])]    
IDdf.insert(2, 'sampleMut', IDdf_sampleMut)

# 2) indel-ID83 Mutation Probabilities
id83ActivitiesPath = '/mnt/e/5_signatureResult34_41samples/42ID83/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/'
IDmutProb = pd.read_csv(id83ActivitiesPath + 'De_Novo_MutationType_Probabilities.txt', sep='\t')   
IDsampleMut = [str(j)+'_' + str(k) for j, k in zip(IDmutProb['Sample Names'], IDmutProb['MutationType'])]
IDmutProb.insert(2, 'sampleMut',IDsampleMut)
IDmutProb = IDmutProb[['sampleMut',	'ID1',	'ID2',	'ID5',	'ID6',	'ID8',	'ID9',	'ID14']]

# 3) merged by the common column sampleMut 
df_IDseqProb = IDdf.merge(IDmutProb)
df_IDseqProb['case_chr_pos_mut'] = df_IDseqProb.apply(lambda row: '_'.join([row['sampleMut'].split('_')[0], str(row['chrs']), str(row['pos']), row['sampleMut'].split('_')[2]]),axis=1)
df_IDseqProb = df_IDseqProb[['samples',	'case_chr_pos_mut', 'chrs',	'pos',	'ID1', 'ID2', 'ID5', 'ID6',	'ID8', 'ID9', 'ID14']]
df_IDseqProb.insert(0, 'caseid', [x.split('_')[0] for x in df_IDseqProb.samples])
df_IDseqProb.insert(0, 'case_chrs_pos', df_IDseqProb['caseid'].astype(str) + '_' + df_IDseqProb['chrs'].astype(str)
                    + '_' + df_IDseqProb['pos'].astype(str))
df_IDseqProb = df_IDseqProb[['case_chrs_pos', 'samples','case_chr_pos_mut',	'ID1',	'ID2',	'ID5',	'ID6',	'ID8',	'ID9',	'ID14']]
df_IDseqProb = df_IDseqProb.replace(sampleOldID, sampleNewID, regex=True)
df_IDseqProb = df_IDseqProb.query('samples != "case12_S12"')
df_IDseqProb = df_IDseqProb.drop_duplicates(subset='case_chr_pos_mut', keep='first')
df_IDseqProb.to_csv('/IPMNPDAC_WGS/Data/Data/sigDPC/step2seqInfoSigOut/41seqInfoIDprob.csv', index=0)
