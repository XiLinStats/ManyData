# install.packages("parallel")
library(parallel)

library(foreach)
library(doSNOW)


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

registerDoSNOW(my.cluster)

# Read in simulated dataset -----------------------------------------------

# data_exp <- readRDS(file.path(data_dir,"/data_exp.rds"))
# data_obs <- readRDS(file.path(data_dir,"/data_obs.rds"))



# Loop --------------------------------------------------------------------

# Number of bootstrap samples
n_boot <- 100

# List of eta
eta_list <- c(seq(0,0.4,0.05),0.5,0.6,0.7,0.8,0.9,1)
# Function to evaluate

set.seed(111)

forms2_obs <- causl:::tidy_formulas(forms_obs[-2], kwd = "cop")
full_form_obs <- causl:::merge_formulas(forms2_obs)
wh_obs <- full_form_obs$wh
# LHS <- lhs(forms2[-length(forms2)])
msks_obs <- ManyData:::masks(forms_obs[-2],family = list(1,1,1),wh_obs)
theta_obs <- causl:::theta(pars = pars_obs, formulas = forms_obs[-2], full_form_obs, kwd = "cop")


forms2_exp <- causl:::tidy_formulas(forms_exp[-2], kwd = "cop")
full_form_exp <- causl:::merge_formulas(forms2_exp)
wh_exp <- full_form_exp$wh
# LHS <- lhs(forms2[-length(forms2)])
msks_exp <- ManyData:::masks(forms_exp[-2],family = list(1,1,1),wh_exp)
theta_exp <- causl:::theta(pars = pars_exp, formulas = forms_exp[-2], full_form_exp, kwd = "cop")



func <- function(){

  data_exp_b <- causalSamp(250,formulas = forms_exp,pars = pars_exp, family = list(1,c(5,5),1,1))
  data_obs_b <- data_obs <- causalSamp(2500,formulas = forms_obs,pars = pars_obs, family = list(1,c(5,5),1,1))


  fit_exp <- fitCausal(dat = data_exp_b,
                     formulas = forms_exp[-2],
                     family = c(1,1,1))

  MLE_exp<- fit_exp$pars$Y$beta
  MLE_exp[c(2,3)] <- fit_exp$pars$Y$beta[c(3,2)]


  result<-matrix(0,1,7)

  for (eta in eta_list){

  theta_sim_raw <- run_MH_MCMC(MLE_exp, 5000, data_obs_b, data_exp_b,
                               msks_obs,theta_obs,
                               msks_exp,theta_exp,MLE_exp, 0.2*diag(4),eta)

  theta_sim_rtnd <- theta_sim_raw[-seq(1:50),]


  lst <- list()
  for (i in 1:nrow(theta_sim_rtnd)) {

    theta2 <- copy(theta_exp)
    theta2[c(3,4,5,6)] <- theta_sim_rtnd[i,]
    msks2 <- copy(msks_exp)
    np <- sum(msks2$beta_m > 0)
    msks2$beta_m[msks_exp$beta_m > 0] <- theta2[seq_len(np)]
    msks2$phi_m[msks_exp$phi_m > 0] <- theta2[-seq_len(np)]
    mm_exp <- model.matrix(full_form_exp$formula, data = data_exp_b)
    ll_i <- ManyData:::llC(data_exp_b[,c("Z","Y")],mm_exp,msks2$beta_m, phi = msks2$phi_m,inCop = c(1,2))
    lst[[i]] <- as.vector(ll_i)
    # print(i)

  }

  mtrx <- do.call(rbind, lst)
  waic_eta <- waic(mtrx)

  result<-rbind(result,c(eta,colMeans(theta_sim_rtnd),mean(data_exp_b$C1),waic_eta$estimates[1]))

  }
  return(result[-1,])

}


start <- Sys.time()

pb <- txtProgressBar(max = n_boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(
  i = 1:n_boot,
  .combine = 'rbind',
  .packages = c('mvtnorm','data.table','coda','ManyData','loo','causl'),
  .options.snow = opts
) %dopar% {
  # .GlobalEnv$forms_obs <- forms_obs
  .GlobalEnv$full_form_obs <- full_form_obs
  # .GlobalEnv$msks_obs <- msks_obs
  # .GlobalEnv$masks <- masks
  # .GlobalEnv$pars_obs <- pars_obs
  # .GlobalEnv$forms_exp <- forms_exp
  .GlobalEnv$full_form_exp <- full_form_exp
  # .GlobalEnv$MLE_exp <- MLE_exp
  x<-func()
  # x<- summary(rnorm(1e7))[3]
  return(x)
}


results <- as.data.table(results)
colnames(results) <- c("eta","beta_0","beta_2","beta_1","beta_3", "meanC1_exp","ELPD_waic")


Sys.time() -start
saveRDS(results,file.path(results_dir,"results_100iter.rds"))
# results<-readRDS(file.path(results_dir,"results_100iter.rds"))

beta_1_true <- 0.8
beta_3_true <- 0.3
# ATE_true <- beta_1_true + beta_3_true * expit(0)

results[,ATE:= beta_1 + beta_3 * meanC1_exp]

summ<-results[,.(.N,beta_1 = mean(beta_1),
                beta_1sd = sd(beta_1),
                beta_3 = mean(beta_3),
                beta_3sd = sd(beta_3),
                ATE = mean(ATE),
                ATE_sd = sd(ATE),
                elpd = mean(ELPD_waic),
                RMSE_beta1 = sqrt(mean((beta_1 - beta_1_true)^2)),
                RMSE_beta3 = sqrt(mean((beta_3 - beta_3_true)^2)),
                RMSE_ATE = sqrt(mean((ATE - (beta_1_true + beta_3_true *meanC1_exp ))^2))) ,eta]



p1 <- ggplot(summ, aes(x = eta)) +
  geom_line(aes(y = -elpd), color = "darkred", size = 1) +
  ggtitle("Estimated ELPD")+
  theme_bw()
p1


p2 <- ggplot(summ, aes(x = eta)) +
  geom_line(aes(y = RMSE_beta1), color = "darkred", size = 1) +
  ggtitle("RMSE beta1")+
  theme_bw()
p2



p3 <- ggplot(summ, aes(x = eta)) +
  geom_line(aes(y = RMSE_beta3), color = "darkred", size = 1) +
  ggtitle("RMSE beta3")+
  theme_bw()
p3


p4 <- ggplot(summ, aes(x = eta)) +
  geom_line(aes(y = RMSE_ATE), color = "darkred", size = 1) +
  ggtitle("RMSE ATE")+
  theme_bw()
p4


grid.arrange(p1,p2,p3,p4,nrow = 2)
