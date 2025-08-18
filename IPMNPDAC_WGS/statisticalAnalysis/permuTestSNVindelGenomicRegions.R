%%R
library("coin")

set.seed(123) #if not, output slight different

# 1) pair-wise nopaired permutation test function
pairwise.permutation.test = 
	function(x, g, data, method = "fdr")
	{n = length(unique(g))
		N = n*(n-1)/2
		d = data.frame(x = x, g = g)
		Z = data.frame(Comparison=rep("A", N),W=rep(NA, N),
									  p.value=rep(NA, N),
									  p.adjust=rep(NA, N),
									  stringsAsFactors=FALSE)
		k=0               
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				k=k+1
				Namea = as.character(unique(g)[i])
				Nameb = as.character(unique(g)[j])
				Datax = subset(d, g==unique(g)[i])
				Datay = subset(d, g==unique(g)[j])
				Dataz = rbind(Datax, Datay)
				Dataz$g2 = factor(Dataz$g)
                #z = independence_test(x ~ g2, data=Dataz, distribution = "asymptotic")
				z = oneway_test(x ~ g2, data=Dataz, distribution = approximate(nresample = 9999))
				P = signif(pvalue(z), digits=4)
				S = signif(statistic(z), digits=4)
				P.adjust = NA                       
				Z[k,] = c(paste0(Namea, "_vs_", Nameb), S, P, P.adjust)
			}
		} 
		Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits = 4) 
		Z
	}

# 2) input dataset for above test
df = read.csv('/IPMNPDAC_WGS/Data/ipmnsnvindel4permutest.csv')
tumorStage <-c('IPMN', 'PDAC')
clNames = names(df)[-1]
combnNumber = dim(combn(length(unique(tumorStage)),2))[2]

# 3) data one step pair-wise test
dfs = data.frame(combnNumber)
for (colnm in clNames){X = df[,colnm]
                       A <-df$tumorStage
                       testDf = data.frame(A,X)
                       w = pairwise.permutation.test(X,A, data = testDf, method='fdr')
                       w = w[,c(1, 3,4)]
                       names(w) <-c(paste0('pairWise', names(w)[1]), 
                                    paste0(colnm,'_',names(w)[2]), 
                                    paste0(colnm,'_',names(w)[3]))
                       dfs <- cbind(dfs, w)
                       }

dfsb = dfs[, !duplicated(colnames(dfs))][,-1]
#write.csv(dfsb, '/mnt/a/famousdata/IPMN_SNV_Indel_pairwiseTestPvalues.csv', row.names=FALSE)
##get pretty print
library(knitr)
library(kableExtra)
library(IRdisplay)
kable(dfsb, "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
  dfsb
