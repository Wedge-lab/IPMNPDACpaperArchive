
import glob
import pandas as pd
import numpy as np
from natsort import natsorted
from natsort import index_natsorted

sampleID = pd.read_csv('/mnt/b/1_3July23VersionDataPlotCode/Figure1/sampleID2025.csv')
sampleNewID = list(sampleID.ID2025)
sampleOldID = list(sampleID.previousID)

## snv-SBS96 seqinfo seq process, e.g. U:AA[C>A]GA -> A[C>A]G
seqinfoPath = '/mnt/e/5_signatureResult34_41samples/41SBS96/output/vcf_files/SNV/'  
outpath = '/IPMNPDAC_WGS/Data/Data/sigDPC/step1seqProcessed/'
sbsDFs = []
for fn in glob.glob(seqinfoPath+'*seqinfo.txt'):
    #combined 22 files for 22 chromsomes
    sbsDf = pd.read_csv(fn, sep='\t', header=None)
    sbsDf.columns = ['samples', 'chrs', 'pos', 'mutationType', 'orientation']
    sbsMutTypeSeq = [a[-8:-1] for a in sbsDf.mutationType]
    sbsDf.insert(4, 'mutationTypeSeq', sbsMutTypeSeq)
    sbsDFs.append(sbsDf)

sampleList = sorted(list(set(sbsDFs[0].samples)))
SBSdf = pd.concat(sbsDFs) 
SBSdf_sampleMut = [str(b)+'_'+str(c) for b, c in zip(SBSdf['samples'], SBSdf['mutationTypeSeq'])]    
SBSdf.insert(2, 'sampleMut', SBSdf_sampleMut)  
SBSdf = SBSdf.replace(sampleOldID, sampleNewID, regex=True)
sampleList = sorted(list(set(SBSdf.samples)))
SBSdf.to_csv(outpath + '41Samples_SBSseqinfo.csv', index=0) #all combined datasets
#########
for samplex in sampleList:
    sampleSBSseq = SBSdf[SBSdf['samples'].str.contains(samplex)]
    sampleSBSseq.insert(0, 'chr_pos',[str(p)+'_'+str(q) for p,q in zip(sampleSBSseq.chrs, sampleSBSseq.pos)])
    sampleSBSseq.to_csv(outpath+ '{}_SBSseqinfo.csv'.format(samplex), index=0)  
#######
