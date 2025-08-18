import sys
import pandas as pd
from natsort import natsort_keygen

def bbDPCtiming(sampleID, ssDPIpath, ssDPCOpath, tmOpath):
    #0) sampleID and path as parameters
    sampleID = sampleID
    ssDPCIfolderpath = ssDPIpath
    dpcoutfolderpath = ssDPCOpath
    outpath = tmOpath

    dpcoutSampleFolderSuffix = '_DPoutput_2000iters_1000burnin_seed123/'

    # 1) raw inputs
    ssdpi_txt = ssDPCIfolderpath + sampleID+'_ssDPI.txt'
    dpcoutCluster_txt = dpcoutfolderpath  + sampleID + dpcoutSampleFolderSuffix + sampleID + '_2000iters_1000burnin_bestClusterInfo.txt'
    dpcout_bed = dpcoutfolderpath + sampleID + dpcoutSampleFolderSuffix + sampleID + '_2000iters_1000burnin_bestConsensusAssignments.bed'

    # 2) working inputs: combine ssDPI and dpcout *bestConsensusAssignments.bed to estimate timing
    ssdpi = pd.read_csv(ssdpi_txt, sep='\t')
    bestCluster = pd.read_csv(dpcoutCluster_txt, sep='\t')
    dpcOut = pd.read_csv(dpcout_bed , sep='\t')

    #ssdpi = ssdpi[['chr','end','subclonal.CN','nMaj1','frac1','nMaj2','frac2','no.chrs.bearing.mut']]
    #version2 add nMin1 and nMin2
    ssdpi = ssdpi[['chr','end','subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2','no.chrs.bearing.mut']]
    ssdpi.insert(0, 'chr_pos', [str(a)+'_'+str(b) for a, b in zip(ssdpi.chr, ssdpi.end)]) 

    dpcOut = dpcOut[['chr','end','cluster']].fillna(0) #fill value 0 as unsp timing
    bestCluster = bestCluster.rename(columns={"cluster.no": "cluster", 'location':'location_CCF'})
    bestCluster = bestCluster[['cluster', 'location_CCF']]
    bestClusterb = pd.DataFrame({'cluster':[0], 'location_CCF':[0]})
    bestClusterNa = pd.concat([bestCluster,bestClusterb])

    dpcOut_CCF = dpcOut.merge(bestClusterNa)
    dpcOut_CCF.insert(0, 'chr_pos', [str(c)+'_'+str(d) for c, d in zip(dpcOut_CCF.chr, dpcOut_CCF.end)])
    dpcOut_CCF = dpcOut_CCF.sort_values(by='chr_pos', key=natsort_keygen())
    dpcOut_CCF = dpcOut_CCF[['chr_pos', 'cluster', 'location_CCF']]

    noCNA_loci = list(set(ssdpi.chr_pos)-set(dpcOut_CCF.chr_pos))
    noCNA_df = pd.DataFrame({'chr_pos':noCNA_loci})
    dpcOut_CCFna = pd.DataFrame(columns=list(dpcOut_CCF)[1:], index=range(noCNA_df.shape[0]))
    dpcOut_CCFNa = pd.concat([noCNA_df, dpcOut_CCFna], axis=1).fillna(-1) #fill value < 0 as unsp timing

    finalDpcOut_CCF = pd.concat([dpcOut_CCF, dpcOut_CCFNa])
    ssdpi_dpcOut = ssdpi.merge(finalDpcOut_CCF)
    ssdpi_dpcOut = ssdpi_dpcOut.rename(columns={'no.chrs.bearing.mut':'noChrsBearingMut'})

    ssdpi_dpcOut.to_csv(outpath+'{}_ssdpi_dpcOutMerge.csv'.format(sampleID), index=0)

    # 3) get dpc_timing
    tms = []
    for i, ccf in enumerate(ssdpi_dpcOut.location_CCF):
        if ccf <= 0: #fill values <=0 as unsp, just for coding selection
            tm = 'unSp'
        elif 0 < ccf < 0.95:
            tm ='subclone'
        elif ssdpi_dpcOut.nMaj1[i] > 1 or ssdpi_dpcOut.nMaj2[i] > 1:
            if ssdpi_dpcOut.noChrsBearingMut[i] > 1:
                tm = 'cloneEarly'
            #version2-add possiblyLate
            elif ssdpi_dpcOut.nMin1[i] == 1 or ssdpi_dpcOut.nMin2[i] == 1 and ssdpi_dpcOut.noChrsBearingMut[i] <= 1:
                tm = "possiblyLate"
            else:
                tm = 'cloneLate'
        else:
            tm = 'cloneNA'
        tms.append(tm)

    ssdpi_dpcOut.insert(ssdpi_dpcOut.shape[1], 'bbDPCtiming', tms)
    ssdpi_dpcOut = ssdpi_dpcOut.sort_values(by='chr_pos', key=natsort_keygen())
    ssdpi_dpcOut.to_csv(outpath+'{}_bbDPCtimg.csv'.format(sampleID), index=0)
    
def doit4me(args):
    bbDPCtiming(args[1],args[2],args[3],args[4])   
    
if __name__ == '__main__':
    doit4me(sys.argv)
    
