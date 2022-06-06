

# Install packages --------------------------------------------------------


# install.packages("parallel")
pkg_list <- c("data.table","survey","purrr","loo","mvtnorm","foreach","doSNOW","coda","doParallel")

for (pkg in pkg_list) {
  if (!(pkg %in% installed.packages())) {
    install.packages(pkg)
  }
}


# install.packages("frugalSim")
# ManyData_path <- file.path("D:/xlin/R_code/Server/code4/ManyData.tar.gz")
# install.packages(ManyData_path, repos = NULL, type = "source")
# devtools::install_github("XiLinStats/ManyData")

library(causl)
library(ManyData)
library(data.table)
library(survey)
library(mvtnorm)
library(doParallel)
library(loo)
library(coda)
library(coda)
library(parallel)
library(foreach)
library(doSNOW)


# Parallelisation ---------------------------------------------------------

parallel::detectCores()
# [1] 12
n.cores <- parallel::detectCores() - 1
# n.cores <- 44
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

# Parameters ------------------------------------------------------------------

bias_1 <- 0.06
bias_2 <- 0.09


forms_obs <- list(Z ~ C1,
                  list(X ~ Z*C1,C1~ 1),
                  Y ~ X*C1,
                  ~ X*C1 #copula
)


# parameters


pars_obs <- list(C1 = list(beta = 0),
                 Z = list(beta = c(0.2,0.6),phi = 1),
                 X = list(beta = c(0.5,0.6,0.1,0.4)),
                 Y = list(beta = c(0.6,0.8 + bias_1 ,0.2,0.3+ bias_2),phi = 1),
                 cop = list(beta = matrix(c(1,0,0,0),nrow = 4)))


forms_exp <- list(Z ~ C1,
                  list(X ~ 1,C1~ 1),
                  Y ~ X*C1,
                  ~ X*C1 #copula
)


# parameters

pars_exp <- list(C1 = list(beta = 0),
                 Z = list(beta = c(0.2,0.6),phi = 1),
                 X = list(beta = c(0)),
                 Y = list(beta = c(0.6,0.8,0.2,0.3),phi = 1),
                 cop = list(beta = matrix(c(1,0,0,0),nrow = 4)))


# Loop --------------------------------------------------------------------


# Number of repeats (number of different datasets)
n_boot <- 100

# List of eta
# eta_list <- c(seq(0,0.2,0.04),seq(0.25,0.4,0.05),seq(0.5,1,0.1))
eta_list <- c(0.2)

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

  # Simulate Data

  data_exp_b <- causalSamp(250,formulas = forms_exp,pars = pars_exp, family = list(1,c(5,5),1,1))
  data_obs_b <- causalSamp(2500,formulas = forms_obs,pars = pars_obs, family = list(1,c(5,5),1,1))


  # Use the MLE fit in the experimental data to initialise the MCMC chain, as well as prior
  fit_exp<-fitCausal(dat = data_exp_b,
                     formulas = forms_exp[-2],
                     family = c(1,1,1)
  )

  MLE_exp <- c(fit_exp$pars$Z$beta,fit_exp$pars$Y$beta[1],fit_exp$pars$Y$beta[3],fit_exp$pars$Y$beta[2],fit_exp$pars$Y$beta[4],
               fit_exp$pars$cop$beta[1], fit_exp$pars$cop$beta[3],fit_exp$pars$cop$beta[2],fit_exp$pars$cop$beta[4],
               fit_exp$pars$Z$phi,fit_exp$pars$Y$phi

  )

  result<-matrix(0,1,9)

  for (eta in eta_list){

    # Get posterior sample using MCMC

    theta_sim_raw <- run_MH_MCMC(MLE_exp, 15000, data_obs_b, data_exp_b,
                                 msks_obs,theta_obs,
                                 msks_exp,theta_exp,MLE_exp, 0.2*diag(12),eta)

    # Burn in
    theta_sim_rtnd <- theta_sim_raw[-seq(1:500),]

    # Calculate Acceptance Ratio
    AR <- mean(c(length(unique(theta_sim_rtnd[,3])),
                 length(unique(theta_sim_rtnd[,4])),
                 length(unique(theta_sim_rtnd[,5])),
                 length(unique(theta_sim_rtnd[,6]))))/nrow(theta_sim_rtnd)


    # Apply thinning
    theta_sim_rtnd <- theta_sim_rtnd[seq(1,nrow(theta_sim_rtnd),5),]

    # Calculate effective sample size
    ESS <- mean(effectiveSize(theta_sim_raw)[c(3,4,5,6)])


    # Calculate ELPD

    lst <- list()
    for (i in 1:nrow(theta_sim_rtnd)) {

      # theta2 <- copy(theta_exp)
      theta2 <- theta_sim_rtnd[i,]
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


    # Output results
    result<-rbind(result,c(eta,AR,ESS, colMeans(theta_sim_rtnd[,c(3,4,5,6)]),mean(data_exp_b$C1),waic_eta$estimates[1]))

  }
  return(result[-1,])

}

# func()


# Run loops in parallel ---------------------------------------------------


start <- Sys.time()

print(paste0("start of parallelisation :",Sys.time()))

results1 <- foreach(
  i = 1:n_boot,
  .combine = 'rbind',
  .packages = c('mvtnorm','data.table','coda','ManyData','loo','causl')

) %dopar% {
  .GlobalEnv$full_form_obs <- full_form_obs
  .GlobalEnv$full_form_exp <- full_form_exp
  x<-func()
  return(x)
}
results1 <- as.data.table(results1)
colnames(results1) <- c("eta","AR","ESS","beta_0","beta_2","beta_1","beta_3", "meanC1_exp","ELPD_waic")
# print(results1)
saveRDS(results1,"data/medium_results.rds")

print(Sys.time()-start)


