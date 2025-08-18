import pandas as pd 

# 1) process Probabilities file by adding column sampleMut: e.g. case10_2_A[C>A]A
sampleID = pd.read_csv('/mnt/b/1_3July23VersionDataPlotCode/Figure1/sampleID2025.csv')
sampleNewID = list(sampleID.ID2025)
sampleOldID = list(sampleID.previousID)

pathsigout = '/mnt/e/5_signatureResult34_41samples/41SBS96/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/'
sigoutseq = pd.read_csv(pathsigout + 'De_Novo_MutationType_Probabilities.txt', sep = '\t')
sigoutseq = sigoutseq.replace(sampleOldID, sampleNewID, regex=True)
sigoutseq.insert(0, 'sampleMut', [str(a) +'_'+str(b) for a,b in zip(sigoutseq.SampleNames, sigoutseq.MutationType)])

# 2) merge seqinfo and SBS probabilities on the column'sampleMut'	
seqInfopath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step1seqProcessed/'
seqInfoSigOut = '/IPMNPDAC_WGS/Data/Data/sigDPC/step2seqInfoSigOut/'
allSampleSeqInfo = pd.read_csv(seqInfopath + '41Samples_SBSseqinfo.csv')
segInfoSBSprob = allSampleSeqInfo.merge(sigoutseq)
colNameSBS = [a for a in list(segInfoSBSprob) if 'SBS' in a]
segInfoSBSprob = segInfoSBSprob[['samples', 'chrs', 'pos', 'sampleMut'] + colNameSBS]
segInfoSBSprob.to_csv(seqInfoSigOut +'41seqInfoSBSprob.csv', index=0)
