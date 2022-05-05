## get param masks
forms2_exp <- causl:::tidy_formulas(forms_exp[-2], kwd = "cop")
full_form_exp <- causl:::merge_formulas(forms2_exp)
wh_exp <- full_form_exp$wh
# LHS <- lhs(forms2[-length(forms2)])
msks_exp <- masks(forms_exp[-2],family = list(1,1,1),wh_exp)

## get model matrix
mm_exp <- model.matrix(full_form_exp$formula, data = data_exp)

theta_exp <- causl:::theta(pars = pars_exp, formulas = forms_exp[-2], full_form_exp, kwd = "cop")


eta_select <- function(eta, n_iter, burn_in, thin){

  theta_sim_raw <- run_MH_MCMC(c(1,0.5), n_iter, dat_obs=data_obs, dat_exp = data_exp, eta = eta)

  accept_rate <- nrow(unique(theta_sim_raw))/nrow(theta_sim_raw)
  theta_sim_rtnd <- theta_sim_raw[-seq(1:burn_in),][seq(1, nrow(theta_sim_raw) - burn_in, by =thin), ]
  ess <- effectiveSize(theta_sim_rtnd)

  MSE_1 <- mean((theta_sim_rtnd[,1] - 0.2)^2)
  MSE_2 <- mean((theta_sim_rtnd[,2] - 0.8)^2)

  lst <- list()


  for (i in 1:nrow(theta_sim_rtnd)) {

    theta2 <- copy(theta_exp)
    theta2[c(6,8)] <- theta_sim_rtnd[i,]

    msks2 <- copy(msks_exp)
    np <- sum(msks2$beta_m > 0)
    msks2$beta_m[msks_exp$beta_m > 0] <- theta2[seq_len(np)]
    msks2$phi_m[msks_exp$phi_m > 0] <- theta2[-seq_len(np)]
    ll_i <- causl:::ll(dat = data_exp[,c(1,5)], mm = mm_exp, beta = msks2$beta_m, phi = msks2$phi_m, inCop = c(1,2),
                       fam_cop = 1, family = list(1,1), link = NULL, par2 = NULL,
                       useC = TRUE)

    lst[[i]] <- ll_i
    printCount(i)

  }

  mtrx <- do.call(rbind, lst)
  loo_eta <- loo(mtrx)
  waic_eta <- waic(mtrx)

  result <- list(loo_eta,waic_eta,MSE_1,MSE_2,ess,accept_rate,mtrx,theta_sim_rtnd)

  return(result)

}

# result<- eta_select(eta = 1, n_iter = 100, burn_in = 20)


eta_list <- seq(0,1,0.1)

elpd_waic <- list()
elpd_waic_se <- list()
elpd_loo <- list()
elpd_loo_se <- list()
MSE_1 <- list()
MSE_2 <- list()
ESS_1 <- list()
ESS_2 <- list()
AR <- list()
mtrx_list <- list()
theta_sim_list <- list()
start_time <- Sys.time()

for (i in 1 : length(eta_list)) {

  eta <- eta_list[[i]]

  print(paste0("The current eta: ",eta))

  result <- eta_select(eta = eta, n_iter = 50000, burn_in = 100, thin = 5)
  elpd_loo <- append(elpd_loo,result[[1]]$estimates[1])
  elpd_loo_se <- append(elpd_loo_se,result[[1]]$estimates[1,2])
  elpd_waic <- append(elpd_waic,result[[2]]$estimates[1])
  elpd_waic_se <- append(elpd_waic_se,result[[2]]$estimates[1,2])
  MSE_1 <- append(MSE_1,result[[3]])
  MSE_2 <- append(MSE_2,result[[4]])
  ESS_1 <- append(ESS_1,result[[5]][1])
  ESS_2 <- append(ESS_2,result[[5]][2])
  AR <- append(AR,result[[6]])
  mtrx_list[[i]] <- result[[7]]
  theta_sim_list[[i]] <- result[[8]]

}

end_time <- Sys.time()
end_time -  start_time


table <- as.data.table(list(eta = eta_list,
                            AR = AR,
                            ESS_1 = ESS_1,
                            ESS_2 = ESS_2,
                            elpd_waic = -unlist(elpd_waic),
                            elpd_loo = -unlist(elpd_loo),
                            elpd_waic_se = unlist(elpd_waic_se),
                            elpd_loo_se = unlist(elpd_loo_se),
                            MSE_1 = unlist(MSE_1),
                            MSE_2 = unlist(MSE_2)))

table[,RMSE_1 := MSE_1^(1/2)]
table[,RMSE_2 := MSE_2^(1/2)]

output_dir <- "D:/xlin/R_code/WIP-SMI/Results/Baseline/Result"
write.csv(apply(table,2,as.character), file.path(output_dir,"table.csv"))
saveRDS(mtrx_list, file.path(output_dir,"mtrx_list.rds"))
saveRDS(theta_sim_list, file.path(output_dir,"theta_sim_list.rds"))

p1 <- ggplot(table, aes(x = eta)) +
  geom_line(aes(y = elpd_loo), color = "darkred", size = 1) +
  ggtitle("Estimated ELPD")+
  theme_bw()

p2 <- ggplot(table, aes(x = eta)) +
  geom_line(aes(y = RMSE_1), color = "darkred", size = 1) +
  ggtitle("RMSE for beta(X)")+
  theme_bw()

p3 <- ggplot(table, aes(x = eta)) +
  geom_line(aes(y = RMSE_2), color = "darkred", size = 1) +
  ggtitle("RMSE for beta(C1)")+
  # geom_line(aes(y = elpd_loo), color ="steelblue", size =2,linetype = "twodash") +
  theme_bw()


grid.arrange(p1, p2, p3, ncol = 1)

loo_compare(loo(mtrx_list[[1]]), loo(mtrx_list[[11]]))
loo_compare(loo(mtrx_list[[1]]), loo(mtrx_list[[2]]))










































