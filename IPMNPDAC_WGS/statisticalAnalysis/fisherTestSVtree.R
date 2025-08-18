
%%R
library("coin")
library(stringr)

set.seed(123)
df = read.csv('/IPMNPDAC_WGS/Data/SV_multiple_single_branches.csv')
X = as.numeric(df$Total)
A = as.factor(df$samples)
z = oneway_test(X~A, data=df,  distribution = approximate(nresample = 9999))
z
