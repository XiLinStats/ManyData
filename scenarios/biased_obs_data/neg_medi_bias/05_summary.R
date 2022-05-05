library(writexl)

bounds_table<- readRDS(file.path(data_dir, "bounds_table.RDS"))

# ATE estimates -----------------------------------------------------------

ATE_fitCausl_exp <- MLE_exp[3] + MLE_exp[4] * mean(data_exp$C1)
ATE_MLE_exp <- MLE_fix_exp[3] + MLE_fix_exp[4] * mean(data_exp$C1)


# bounds_table[,ATE_exp_l := beta_1l + beta_3l * mean(data_exp$C1)]
# bounds_table[,ATE_exp_u := beta_1u + beta_3u * mean(data_exp$C1)]

bounds_table[,ATE_fitCausl_exp_ind := 1*(ATE_fitCausl_exp < ATE_exp_u & ATE_fitCausl_exp > ATE_exp_l)]
bounds_table[,.(cover_prob = mean(ATE_fitCausl_exp_ind)),eta]

bounds_table[,ATE_MLE_exp_ind := 1*(ATE_MLE_exp < ATE_exp_u & ATE_MLE_exp > ATE_exp_l)]
bounds_table[,.(cover_prob = mean(ATE_MLE_exp_ind)),eta]

summ<-bounds_table[,.(AR = mean(AR),
                ESS = mean(ESS),
                cover_prob_fitCausal = mean(ATE_fitCausl_exp_ind),
                cover_prob_fix = mean(ATE_MLE_exp_ind)),eta]

write_xlsx(list(R_output = summ),file.path(results_dir,"summary_large_bias.xlsx"))

