import pandas as pd
from scipy.stats import ttest_ind
datapath = '/mnt/b/1_3July23VersionDataPlotCode/NC_revision/excelSubmited/'
purityData = pd.read_csv(datapath + '573583_0_data_set_10134107_spmdqf.csv')
purityData = purityData[['sampleName', 'cellularity']]

multipleBranch = ['case10','case3','case16','case2','case9','case7','case6']
singleBranch = ['case4','case12','case13', 'case15']
sampleNames = list(purityData.sampleName)

multipleBranchSample = [x for x in sampleNames if x.split('_')[0] in multipleBranch]
singleBranchSample = [x for x in sampleNames if x.split('_')[0] in singleBranch]

multipleBranchPurity = purityData.query('sampleName==@multipleBranchSample')
singleBranchPurity = purityData.query('sampleName==@singleBranchSample')
mean_multipleBranchPurity = multipleBranchPurity['cellularity'].mean()
mean_singleBranchPurity = singleBranchPurity['cellularity'].mean()
multipleBranchPurityV = list(multipleBranchPurity.cellularity)
singleBranchPurityV = list(singleBranchPurity.cellularity)

t_stat, p_value = ttest_ind(multipleBranchPurityV , singleBranchPurityV)
