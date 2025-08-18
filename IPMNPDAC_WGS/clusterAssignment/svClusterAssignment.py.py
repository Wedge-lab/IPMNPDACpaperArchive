import os
from glob import glob
import pandas as pd
import pyranges as pr

# a) function for SV matching cluster
def msClusterSV(msDPCdata, svdata, pathout):
    # 1) Load data
    caseMsDPC = pd.read_csv(msDPCdata)
    caseID = msDPCdata.split('/')[-1].split('_')[0]
    caseMsDPC['chrom'] = 'chr' + caseMsDPC['chr'].astype(str)
    caseMsDPC = caseMsDPC.rename(columns={"most.likely.cluster": "msCluster"})

    allsv = pd.read_csv(svdata)
    caseSV = allsv[allsv['samples'].str.contains(caseID)]

    # 2) Convert to PyRanges
     # SNV is a single base interval
    snv = pr.PyRanges(pd.DataFrame({"Chromosome": caseMsDPC['chrom'],"Start": caseMsDPC['pos'],"End": caseMsDPC['pos']+1,  
            "msCluster": caseMsDPC['msCluster']}))

    sv = pr.PyRanges(pd.DataFrame({"Chromosome": caseSV['CHROM'],"Start": caseSV['START'],"End": caseSV['END'],
                                   "typeTumorSample": caseSV['typeTumorSample'],"samples": caseSV['samples'],"SVLEN": caseSV['SVLEN'],
                                   "svclass": caseSV['svclass'],"CHR2": caseSV['CHR2'],"rgn": caseSV['rgn'],"gene": caseSV['gene'],
                                   "tid": caseSV['tid'],"PASS": caseSV['PASS']}))

    #3)  Perform overlap join (very fast C implementation)
    joined = snv.join(sv)
    # Convert back to pandas
    result = joined.df.rename(columns={"Start": "snvPos"})
    result.to_csv(f"{pathout}{caseID}_snvSVCluster.csv", index=False)

# b)Batch execution
def allFolderFileRun():
    pathinput = os.path.expanduser('/IPMNPDAC_WGS/Data/svDriverCluster/CasesCHROMposClusts/')
    sVdatapath = os.path.expanduser('/IPMNPDAC_WGS/Data/svDriverCluster/all41SVs.csv')
    outpathx = os.path.expanduser('/IPMNPDAC_WGS/Data/svDriverCluster/')
    for filex in glob(pathinput + '*_chrPosCluster.csv'):   
        msClusterSV(filex, sVdatapath, outpathx)

if __name__=="__main__":
    allFolderFileRun() 
