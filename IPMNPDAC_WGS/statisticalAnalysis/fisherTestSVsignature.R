
library(coin)

data <- data.frame(
  tumourSample = factor(rep(c("IPMN", "PDAC"), times = c(22, 19))),
  SV2 = factor(c(rep(c("positive", "negative"), times = c(22, 0)),
                      rep(c("positive", "negative"), times = c(17, 2)))),
  SV4 = factor(c(rep(c("positive", "negative"), times = c(0, 22)),
                      rep(c("positive", "negative"), times = c(7, 12)))),
  SV5 = factor(c(rep(c("positive", "negative"), times = c(5, 17)),
                      rep(c("positive", "negative"), times = c(5, 14)))),
  SV7 = factor(c(rep(c("positive", "negative"), times = c(4, 18)),
                      rep(c("positive", "negative"), times = c(9, 10)))),  
  SV9 = factor(c(rep(c("positive", "negative"), times = c(1, 21)),
                       rep(c("positive", "negative"), times = c(5, 14)))) 
)
# save data back
#write.csv(data,'/outpath/SVSigfisherTest.csv', row.names=FALSE)

# Function to run Fisher's Exact Test for each outcome and extract p-value
run_fisher_test <- function(outcome) {
  test_result <- independence_test(reformulate("tumourSample", outcome), data = data, distribution = "exact")
  p_value <- pvalue(test_result)
  return(p_value)
}

outcome_columns <- names(data)[-1]
p_values <- sapply(outcome_columns, run_fisher_test)
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame with the results
results_df <- data.frame(
  Outcome = outcome_columns,
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

print(results_df)
