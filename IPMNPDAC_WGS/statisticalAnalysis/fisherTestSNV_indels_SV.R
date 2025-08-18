
%%R
library(coin)

set.seed(1234)

snv_indel_sv = read.csv('/IPMNPDAC_WGS/Data/snvindelsvTotalnb4StatTest.csv')

snv_indel_sv = snv_indel_sv[,c('tumourType', 'total_indels', 'total_SNVs', 'total_SV')]
#X = as.numeric(tmbdf$snv_indel_TMB)
A = as.factor(snv_indel_sv$tumourType)
Xsnv = as.numeric(snv_indel_sv$total_SNVs)
Xindel = as.numeric(snv_indel_sv$total_indels)
Xsv = as.numeric(snv_indel_sv$total_SV)

snvtest = oneway_test(Xsnv~A, data=snv_indel_sv,  distribution = approximate(nresample = 9999))   
indeltest = oneway_test(Xindel~A, data=snv_indel_sv,  distribution = approximate(nresample = 9999))
svtest = oneway_test(Xsv~A, data=snv_indel_sv,  distribution = approximate(nresample = 9999))
print(snvtest)
print(indeltest)
print(svtest)
