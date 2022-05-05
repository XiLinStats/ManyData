# install.packages("parallel")
library(parallel)

library(foreach)

# Parallelisation ---------------------------------------------------------

parallel::detectCores()
# [1] 12
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()



# Start the loop ----------------------------------------------------------

# True parameters
# theta_true <- c(1.0485330)

# Number of bootstrap samples
n_boot <- 250

# List of eta
eta_list <- c(seq(0,0.4,0.04),0.5,0.6,0.7,0.8,0.9,1)
# Function to evaluate
# eta_list <- c(0)

set.seed(111)

forms2_obs <- causl:::tidy_formulas(forms_obs[-2], kwd = "cop")
full_form_obs <- causl:::merge_formulas(forms2_obs)
wh_obs <- full_form_obs$wh
# LHS <- lhs(forms2[-length(forms2)])
msks_obs <- masks(forms_obs[-2],family = list(1,1,1),wh_obs)
theta_obs <- causl:::theta(pars = pars_obs, formulas = forms_obs[-2], full_form_obs, kwd = "cop")


forms2_exp <- causl:::tidy_formulas(forms_exp[-2], kwd = "cop")
full_form_exp <- causl:::merge_formulas(forms2_exp)
wh_exp <- full_form_exp$wh
# LHS <- lhs(forms2[-length(forms2)])
msks_exp <- masks(forms_exp[-2],family = list(1,1,1),wh_exp)
theta_exp <- causl:::theta(pars = pars_exp, formulas = forms_exp[-2], full_form_exp, kwd = "cop")

data_exp_mean_C1 <- mean(data_exp$C1)
data_obs_mean_C1 <- mean(data_obs$C1)
data_comb_mean_C1 <- mean(data_comb$C1)

func <- function(eta) {

  bound <- matrix(0,1,13)
  data_exp_b <- data_exp[sample(nrow(data_exp), nrow(data_exp),replace = T), ]
  data_obs_b <- data_obs[sample(nrow(data_obs), nrow(data_obs),replace = T), ]

  theta_sim_raw <- run_MH_MCMC(MLE_exp, 5000, data_obs_b, data_exp_b,
                                          msks_obs,theta_obs,
                                          msks_exp,theta_exp,MLE_exp, 0.2*diag(4),eta)

  theta_sim_rtnd <- theta_sim_raw[-seq(1:50),]
  #
  beta_1 <- quantile(theta_sim_rtnd[,3],probs = c(0.02,0.975))
  beta_3 <- quantile(theta_sim_rtnd[,4],probs = c(0.02,0.975))
  ATE_exp <- quantile(theta_sim_rtnd[,3] + theta_sim_rtnd[,4] * data_exp_mean_C1,probs = c(0.02,0.975))
  ATE_obs <- quantile(theta_sim_rtnd[,3] + theta_sim_rtnd[,4] * data_obs_mean_C1,probs = c(0.02,0.975))
  ATE_comb <- quantile(theta_sim_rtnd[,3] + theta_sim_rtnd[,4] * data_comb_mean_C1,probs = c(0.02,0.975))

  ESS <- mean(effectiveSize(theta_sim_raw)[c(3,4)])
  AR <- length(unique(theta_sim_rtnd))/length(theta_sim_rtnd)

  bound[1,] <- c(eta,ESS,AR,beta_1,beta_3,ATE_exp,ATE_obs,ATE_comb)

  return(bound)
}
# Initialise
# cover_prob <- matrix(0,length(eta_list),2)
start_time <- Sys.time()

# Loop


bounds_table <- matrix(0,1,13)
for (i in 1 : length(eta_list)) {

  eta <- eta_list[[i]]
  print(paste0("Current Eta:",eta))

  final <- foreach(
    i = 1:n_boot,
    .combine = 'rbind',
    .packages = c('mvtnorm','data.table','coda','ManyData')
  ) %dopar% {
    .GlobalEnv$forms_obs <- forms_obs
    .GlobalEnv$full_form_obs <- full_form_obs
    .GlobalEnv$msks_obs <- msks_obs
    .GlobalEnv$masks <- masks
    .GlobalEnv$pars_obs <- pars_obs
    .GlobalEnv$forms_exp <- forms_exp
    .GlobalEnv$full_form_exp <- full_form_exp
    .GlobalEnv$pars_exp <- pars_exp
    x<-func(eta)
  }
  # cover_prob[i,] <- c(eta,colMeans(final))
  bounds_table <- rbind(bounds_table,final)
}


bounds_table <- as.data.table(bounds_table)
colnames(bounds_table) <- c("eta","ESS","AR","beta_1l","beta_1u","beta_3l","beta_3u",
                            "ATE_exp_l","ATE_exp_u",
                            "ATE_obs_l","ATE_obs_u",
                            "ATE_comb_l","ATE_comb_u")
Sys.time() -  start_time

saveRDS(bounds_table[-1],file.path(data_dir,"bounds_table.RDS"))


# bounds_table <- rbind(bounds_table,bounds_table_orig[eta!= 0])


























