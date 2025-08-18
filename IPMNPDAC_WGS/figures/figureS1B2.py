import pandas as pd
from natsort import natsort_keygen
import matplotlib.pyplot as plt

from matplotlib import rcParams
plt.rcParams['figure.dpi'] = 300
rcParams['font.family'] = 'Arial'

font = {'family':'monospace','weight':'normal', 'size':12}
#plt.rc('font', **font)

mutpath = "/IPMNPDAC_WGS/Data/"
krasCCF_VAF = pd.read_csv(mutpath + 'kras_vaf_ccf.csv')
krasCCF_VAF = krasCCF_VAF.sort_values(by='typeTumosample', key = natsort_keygen())
krasCCF_VAF['typeTumosample'] = [a[3:] for a in krasCCF_VAF.typeTumosample]

fig, axes = plt.subplots(1,1, figsize=(12,4), sharex=True)
axes.plot(krasCCF_VAF.typeTumosample, krasCCF_VAF.location_CCF,'-o', c='blue',  mfc='black', label='CCF')
axes.plot(krasCCF_VAF.typeTumosample,krasCCF_VAF.cellularity, '-o',c='red', mfc='black', label='Cellularity')
axes.plot(krasCCF_VAF.typeTumosample,krasCCF_VAF.VAF, '-o',c='green', mfc='black', label='VAF') 
axes.set_ylabel('Values')
#axes.set_xtickslabel('Values')
plt.setp(axes.get_xticklabels(), rotation=30, ha="right", fontsize=10)
plt.legend(loc="center left")
plt.show()
