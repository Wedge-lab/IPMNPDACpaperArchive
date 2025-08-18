%%R
library(coin)

tmbdf = read.csv('/IPMNPDAC_WGS/Data/tmbplotPermutest.csv')
X = as.numeric(tmbdf$snv_indel_TMB)
A = as.factor(tmbdf$tumorStage)
tmbtest = oneway_test(X~A, data=tmbdf,  distribution = approximate(nresample = 9999))
tmbtest             
