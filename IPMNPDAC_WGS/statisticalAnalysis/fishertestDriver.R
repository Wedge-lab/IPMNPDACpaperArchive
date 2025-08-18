%%R

library(coin)

data <- data.frame(
  tumourSample = factor(rep(c("IPMN", "PDAC"), times = c(22, 19))),
  KRAS_snvIndel = factor(c(rep(c("positive", "negative"), times = c(10, 12)),
                      rep(c("positive", "negative"), times = c(12, 7)))),
  GNAS_snvIndel = factor(c(rep(c("positive", "negative"), times = c(8, 14)),
                      rep(c("positive", "negative"), times = c(1, 18)))),
  LRP1B_snvIndel = factor(c(rep(c("positive", "negative"), times = c(0, 22)),
                      rep(c("positive", "negative"), times = c(7, 12)))),
  TP53_snvIndel = factor(c(rep(c("positive", "negative"), times = c(2, 20)),
                      rep(c("positive", "negative"), times = c(13, 6)))),  
  U2AF1_los = factor(c(rep(c("positive", "negative"), times = c(4, 18)),
                       rep(c("positive", "negative"), times = c(13, 6)))),
  RNF43_los = factor(c(rep(c("positive", "negative"), times = c(4, 18)),
                       rep(c("positive", "negative"), times = c(12, 7)))),
  los_12q21_31 = factor(c(rep(c("positive", "negative"), times = c(1, 21)),
                       rep(c("positive", "negative"), times = c(10, 9)))),
  los_21q22_3 = factor(c(rep(c("positive", "negative"), times = c(4, 18)),
                       rep(c("positive", "negative"), times = c(13, 6)))),
  los_17q22 = factor(c(rep(c("positive", "negative"), times = c(4, 18)),
                       rep(c("positive", "negative"), times = c(12, 7)))),
  los_6q27 = factor(c(rep(c("positive", "negative"), times = c(8, 14)),
                       rep(c("positive", "negative"), times = c(15, 4)))),
  gain_8q24_3 = factor(c(rep(c("positive", "negative"), times = c(6, 16)),
                       rep(c("positive", "negative"), times = c(16, 3)))),
  gain_14q32 = factor(c(rep(c("positive", "negative"), times = c(15, 7)),
                       rep(c("positive", "negative"), times = c(19, 0)))), 
  gain_17p11 = factor(c(rep(c("positive", "negative"), times = c(19,3)),
                       rep(c("positive", "negative"), times = c(9, 10)))),
  gain_8q23 = factor(c(rep(c("positive", "negative"), times = c(6, 16)),
                       rep(c("positive", "negative"), times = c(12, 7)))),
  gain_11p11 = factor(c(rep(c("positive", "negative"), times = c(3, 19)),
                       rep(c("positive", "negative"), times = c(9, 10))))

)
# save data back
write.csv(data,'/IPMNPDAC_WGS/Data/allIPMNdata4fisherTest.csv', row.names=FALSE)

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

#write.csv(results_df,'/outpath/allIPMNdata4fisherTestP_values.csv', row.names=FALSE)
print(results_df)
