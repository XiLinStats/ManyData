# load packages
library(causl)
library(ManyData)
library(data.table)
library(survey)
library(mvtnorm)
library(doParallel)
library(loo)
library(coda)
library(parallel)
library(foreach)
library(MASS)


# Parallelization ---------------------------------------------------------

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


# Helper functions --------------------------------------------------------

create_strata <- function(K,data_exp_b, data_obs_b){
  q <- seq(0,1,2/K)[-1]
  table <- as.data.table(rbind(as.data.table(list(C5 = 0,
                                                  bin = q,
                                                  C1 = quantile(data_obs_b[C5 == 0]$C1,q))),
                               as.data.table(list(C5 = 1,
                                                  bin = q,
                                                  C1 = quantile(data_obs_b[C5 == 1]$C1,q)))))
  table[bin == 1, C1 := Inf][,K := 1:.N]

  # stratify data
  setkey(table,C5,C1)
  setkey(data_exp_b,C5,C1)
  setkey(data_obs_b,C5,C1)

  data_exp_b <- table[data_exp_b, roll = -Inf]
  data_obs_b <- table[data_obs_b, roll = -Inf]

  return(list(data_exp_b = data_exp_b, data_obs_b = data_obs_b))

}


calc_treat_eff <- function(beta,theta,dat,formula,vars){
  # Create mode matrices
  dat_1 <- copy(dat)
  dat_1[,X := 1]

  dat_0 <- copy(dat)
  dat_0[,X := 0]

  # Create beta matrix
  np <- sum(beta > 0)
  beta[beta > 0] <- theta[seq_len(np)]
  mm_1 <- model.matrix(formula, data = dat_1)
  mm_0 <- model.matrix(formula, data = dat_0)

  # Predict
  eta_1 <- mm_1 %*% beta
  Y1 <- eta_1[,match("Y",vars)]


  eta_0 <- mm_0 %*% beta
  Y0 <- eta_0[,match("Y",vars)]

  eff <- Y1 - Y0
  return(eff)
}

calc_ATE_variance <- function(beta,fit, dat){
  np <- sum(beta > 0)
  var <- 0
  for (j in 1:nrow(dat)) {
    ATE_vec <- as.vector(rep(0,np))
    ATE_vec[11:16] <- c(1,as.numeric(dat[j,c("C1","C2","C3","C4","C5")]))
    var <- var + t(ATE_vec) %*% fit$sandwich[1:np,1:np] %*% ATE_vec
  }
  return(var/nrow(dat))
}


# Parameters --------------------------------------------------------------

n_o <- 2500
n_e <- 250

# Number of iterations (different data)
n_boot <- 500

# List of eta
eta_list <- c(seq(0,0.2,0.04),seq(0.25,1,0.05))

# List of bias
bias_list <- c(seq(0,0.4,0.025),seq(0.45,0.5,0.05),seq(0.6,1,0.1))

# Number of stata used for the estimator in Rosenman et al. (2020)
K <- 10

# Causal model parameterization
family <- list(rep(1,2), c(5,1,1,1,5,5), 1, 1)
forms_obs <- forms_exp <- list(c(Z1 ~ C1,Z2 ~ C4),
                               list(X ~ Z1*C1 + C5, C1 ~ 1, C2 ~ 1, C3 ~ 1, C4 ~ 1, C5 ~ 1),
                               Y ~ C1+ C2+ C3 +C4+ C5+ X + X:C1+ X:C2+ X:C3+ X:C4+ X:C5,
                               ~ X)
forms_exp[[2]][1] <- list(X ~ 1)

pars_exp <- list(Z1 = list(beta = c(0,1),phi = 1),
                 Z2 = list(beta = c(0,1),phi = 1),
                 C1 = list(beta = 0,phi = 1),
                 C2 = list(beta = 0,phi = 1),
                 C3 = list(beta = 0,phi = 1),
                 C4 = list(beta = 0),
                 C5 = list(beta = 0),
                 X = list(beta = c(0.5,0.6,0.1,0.4,0.5)),
                 Y = list(beta = seq(0.1,1.2,0.1),phi = 1),
                 cop = list(beta = matrix(c(1,1,
                                            1,1,
                                            1,1),nrow = 2)))
pars_obs <- pars_exp
pars_exp$X <- list(beta = 0)

forms2 <- list(obs = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"),
               exp = causl:::tidy_formulas(unlist(forms_obs[-2]), kwd = "cop"))


full_form <- list(obs =  causl:::merge_formulas(forms2$obs),
                  exp =  causl:::merge_formulas(forms2$exp))

msks <- list(obs = causl:::masks(forms2$obs,family = c(rep(1,6),1),full_form$obs$wh),
             exp = causl:::masks(forms2$exp,family = c(rep(1,6),1),full_form$exp$wh))

theta <- list(obs = causl:::get_theta(pars = pars_obs, formulas = unlist(forms2$obs), full_form$obs, kwd = "cop"),
              exp = causl:::get_theta(pars = pars_exp, formulas = unlist(forms2$exp), full_form$exp, kwd = "cop"))
vars <- causl:::lhs(unlist(forms2$obs[-length(forms2$obs)]))


# Simulation code ---------------------------------------------------------

set.seed(123)

run_sims <- function(include.oberst = TRUE, include.rosenman = TRUE){

  res <- matrix(0,1,20)

  for (i in 1:length(bias_list)) {

    tryCatch({
      bias <- bias_list[i]

      # generate data
      pars_obs$Y$beta[7] <- pars_exp$Y$beta[7] + bias

      positivity_ind <- FALSE

      while (!positivity_ind) {

        data_obs_b <- as.data.table(causl:::rfrugalParam(n = n_o, formulas = forms_obs, family = family, pars = pars_obs))
        data_exp_b <- as.data.table(causl:::rfrugalParam(n = n_e, formulas = forms_exp, family = family, pars = pars_exp))

        # Stratify

        strata_out <- create_strata(K,data_exp_b,data_obs_b)
        data_exp_b <- strata_out$data_exp_b
        data_obs_b <- strata_out$data_obs_b

        positivity_ind <- !(anyNA(dcast(data_obs_b[,.N,.(X,K)],K ~ paste0("Y_",X)))) & !(anyNA(dcast(data_exp_b[,.N,.(X,K)],K ~ paste0("Y_",X))))

      }


      # Get estimates from both data separately
      fit_exp <- fitCausal(dat = data_exp_b, formulas = unlist(forms2$exp), family = unlist(family[-2]))

      fit_obs <- fitCausal(dat = data_obs_b, formulas = unlist(forms2$obs), family = unlist(family[-2]))

      mm <- list(obs = fit_obs$mm, exp = fit_exp$mm)


      if (include.oberst == TRUE) {
        # Estimator in Oberst et al. (2022)
        np <- sum(msks$exp$beta_m > 0)
        theta_e <- mean(calc_treat_eff(beta = msks$exp$beta_m,theta = fit_exp$par,dat = data_exp_b,formula =  full_form$exp$formula,vars = vars))
        theta_o <- mean(calc_treat_eff(beta = msks$exp$beta_m,theta = fit_obs$par,dat = data_exp_b,formula =  full_form$exp$formula,vars = vars))
        var_e <- calc_ATE_variance(msks$exp$beta_m,fit_exp, data_exp_b)
        var_o <-  calc_ATE_variance(msks$exp$beta_m,fit_obs, data_exp_b)

        ATE_oberst <- obserst(theta_e, theta_o, var_e, var_o)

      } else {
        ATE_oberst <- list(theta_oberst = NA, lambda = NA)
      }

      # Our proposed estimator

      elpd <- -Inf
      opt_samples <- NA
      opt_eta <- NA

      for (eta in eta_list) {
        tryCatch(
          {samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta, n_sample = 2000)
          elpd_sample <- calculate_elpd(samples,data_exp_b[,vars,with = F],mm$exp, msks = msks$exp,method = "WAIC", inCop = 1:3)

          if (elpd < elpd_sample$elpd_waic) {
            elpd <- elpd_sample$elpd_waic
            opt_samples <- samples
            opt_eta <- eta

          }
          }
          ,error = function(e) {cat("ERROR : Bias = ",bias, "eta = ", eta, " skipped", "\n")})
      }

      # There is a chance that the combined covariance is not invertible, theorectically we can then switch back to MCMC but
      # in the simulation, we will just skip to the next bias value

      if (is.na(opt_samples)) {
        next
      }else{
        # Get strata estimates from exp and obs data
        ps_o <- predict(glm(X ~ Z1*C1+C5, family = binomial, data = data_obs_b), type = "response")
        data_obs_b[, ':='(ps = ps_o,
                          weight = X/ps_o + (1 - X)/(1 - ps_o))]

        ps_e <- predict(glm(X ~ Z1*C1+C5, family = binomial, data = data_exp_b), type = "response")
        data_exp_b[, ':='(ps = ps_e,
                          weight = X/ps_e + (1 - X)/(1 - ps_e))]

        # tau_k_e2 <- dcast(data_exp_b[,.(tau = mean(Y)),.(X,K)],K ~ paste0("Y_",X),value.var = "tau")[,tau := Y_1 - Y_0]
        tau_k_e <- dcast(data_exp_b[,.(tau = weighted.mean(Y,w = weight)),.(X,K)],K ~ paste0("Y_",X),value.var = "tau")[,tau := Y_1 - Y_0]
        tau_k_o <- dcast(data_obs_b[,.(tau = weighted.mean(Y,w = weight)),.(X,K)],K ~ paste0("Y_",X),value.var = "tau")[,tau := Y_1 - Y_0]

        tau_e <- as.vector(tau_k_e$tau)
        tau_o <- as.vector(tau_k_o$tau)

        # calculate variance matrix based on the experimental data
        sigma2_k_e <- dcast(data_exp_b[,.(sig = mean((Y - mean(Y))^2)/.N),.(X,K)],K ~ paste0("sigma_",X))[,sigma_sq := sigma_1 + sigma_0]
        sigma2_e <- as.vector(sigma2_k_e$sigma_sq) * diag(K)

        # estimators in Strawderman (1991)
        gsar1 <- tryCatch(gsar1(tau_e, tau_o, Sigma_e = sigma2_e),error = function(e) {return(NA)})
        gsar2 <- tryCatch(gsar2(tau_e, tau_o, Sigma_e = sigma2_e),error = function(e) {return(NA)})

        # estimators in Rosenman et al. (2020)
        kappa1 <- rbob_kappa1(tau_e, tau_o, sigma2_e = sigma2_k_e$sigma_sq)
        kappa2 <- rbob_kappa2(tau_e, tau_o, Sigma_e = sigma2_e)

        # calculate the CATE estimators (true, ours and just using data_exp)
        data_exp_b$eff_e <- calc_treat_eff(beta = msks$exp$beta_m,theta = fit_exp$par,dat = data_exp_b,formula =  full_form$exp$formula,vars = vars)
        data_exp_b$eff_true <- calc_treat_eff(beta = msks$exp$beta_m,theta = theta$exp,dat = data_exp_b,formula =  full_form$exp$formula,vars = vars)
        data_exp_b$eff_eta <- calc_treat_eff(beta = msks$exp$beta_m,theta = colMeans(opt_samples),dat = data_exp_b,formula =  full_form$exp$formula,vars = vars)

        ATE_true <- mean(data_exp_b$eff_true)
        theta_e <- mean(data_exp_b$eff_e)
        ATE_eta <- mean(data_exp_b$eff_eta)

        # Implement method in Kallus et al. (2018)
        mm3 <- as.matrix(cbind(rep(1, n_e),data_exp_b[,c("C1","C2","C3","C4","C5","Z1","Z2")]))
        lm_obs_t <- glm(Y ~ C1+ C2+ C3 +C4+ C5 + Z1 + Z2, family = gaussian, data = data_obs_b[X == 1])
        lm_obs_c <- glm(Y ~ C1+ C2+ C3 +C4+ C5 + Z1 + Z2, family = gaussian, data = data_obs_b[X == 0])
        w_coeff <- lm_obs_t$coefficients - lm_obs_c$coefficients

        data_exp_b[, w_hat := mm3 %*% as.vector(w_coeff)]

        # ps_e <- predict(glm(X ~ Z1*C1, family = binomial, data = data_exp_b), type = "response")
        data_exp_b[, ':='(ps_e = ps_e,q = X/ps_e - (1 - X)/(1 - ps_e))]
        data_exp_b[, resid := q * Y - w_hat]

        eta_model <- glm(resid ~  C1+ C2+ C3 +C4+ C5 + Z1 + Z2, family = gaussian, data = data_exp_b)
        data_exp_b[,eff_kallus := w_hat + eta_model$fitted.values]

        # average unit-level causal effects to get the strata-level effect
        summ <- data_exp_b[,.(CATE_true = mean(eff_true),
                              CATE_eta = mean(eff_eta),
                              CATE_e = mean(eff_e),
                              CATE_kallus = mean(eff_kallus)
        ),K]

        # calculate MSE (CATE)
        CATE_SE_gsar1 <- sum((summ$CATE_true - gsar1)^2)
        CATE_SE_gsar2 <- sum((summ$CATE_true - gsar2)^2)
        CATE_SE_kappa1 <- sum((summ$CATE_true - kappa1)^2)
        CATE_SE_kappa2 <- sum((summ$CATE_true - kappa2)^2)
        CATE_SE_eta <- sum((summ$CATE_true - summ$CATE_eta)^2)
        CATE_SE_kallus <- sum((summ$CATE_true - summ$CATE_kallus)^2)

        CATE_SE_e <- sum((summ$CATE_true - tau_e)^2)
        CATE_SE_e_param <- sum((summ$CATE_true - summ$CATE_e)^2)


        # reweight for ATE
        wt <- data_exp_b[,.N,K]$N
        ATE_gsar1 <- tryCatch(weighted.mean(gsar1,wt),error = function(e) {return(NA)})
        ATE_gsar2 <- tryCatch(weighted.mean(gsar2,wt),error = function(e) {return(NA)})
        ATE_kappa1 <- tryCatch(weighted.mean(kappa1,wt),error = function(e) {return(NA)})
        ATE_kappa2 <- tryCatch(weighted.mean(kappa2,wt),error = function(e) {return(NA)})
        ATE_kallus <- mean(data_exp_b$eff_kallus)

        res <- rbind(res, c(bias, ATE_true, theta_e, ATE_oberst$lambda,ATE_oberst$theta_oberst, opt_eta, ATE_eta,
                            ATE_gsar1,ATE_gsar2,ATE_kappa1,ATE_kappa2, ATE_kallus,
                            CATE_SE_gsar1, CATE_SE_gsar2, CATE_SE_kappa1, CATE_SE_kappa2, CATE_SE_eta, CATE_SE_e,CATE_SE_kallus,
                            CATE_SE_e_param))
      }
    }, error = function(e) {cat("ERROR : Bias = ",bias, " skipped", "\n")})

  }
  return(res[-1,])

}

# Run simulations in parallel ---------------------------------------------

start <- Sys.time()

print(paste0("start of parallelisation :",Sys.time()))

results1 <- foreach(
  i = 1:n_boot,
  .combine = 'rbind',
  .packages = c('mvtnorm','data.table','coda','ManyData','loo','causl','MASS','doParallel')
) %dopar% {
  .GlobalEnv$gsar1 <- gsar1
  .GlobalEnv$gsar2 <- gsar2
  .GlobalEnv$rbob_kappa1 <- rbob_kappa1
  .GlobalEnv$rbob_kappa2 <- rbob_kappa2
  .GlobalEnv$create_strata <- create_strata
  .GlobalEnv$calc_treat_eff <- calc_treat_eff
  .GlobalEnv$calc_ATE_variance <- calc_ATE_variance
  .GlobalEnv$family <- family
  .GlobalEnv$forms_obs <- forms_obs
  .GlobalEnv$forms_exp <- forms_exp
  .GlobalEnv$pars_exp <- pars_exp
  .GlobalEnv$pars_obs <- pars_obs
  .GlobalEnv$forms2 <- forms2
  .GlobalEnv$full_form <- full_form
  .GlobalEnv$msks <- msks
  .GlobalEnv$theta <- theta
  .GlobalEnv$vars <- vars
  x <- run_sims()
  return(x)
}
print(Sys.time() - start)

results1 <- as.data.table(results1)
colnames(results1) <- c("bias", "ATE_true","theta_e", "oberst_lambda","ATE_oberst","opt_eta","ATE_eta",
                        "ATE_gsar1","ATE_gsar2","ATE_kappa1","ATE_kappa2","ATE_kallus",
                        "CATE_SE_gsar1", "CATE_SE_gsar2", "CATE_SE_kappa1", "CATE_SE_kappa2", "CATE_SE_eta", "CATE_SE_e","CATE_SE_kallus","CATE_SE_e_param")
saveRDS(results1,"data/sim_v10.rds")




# Plotting ----------------------------------------------------------------

library(ggplot2)
library(data.table)

dat <- results1

summ2 <- dat[,.(.N,CATE_MSE_gsar1 = mean(CATE_SE_gsar1),
                CATE_MSE_gsar2 = mean(CATE_SE_gsar2 ),
                CATE_MSE_kappa1 = mean(CATE_SE_kappa1 ),
                CATE_MSE_kappa2 = mean(CATE_SE_kappa2 ),
                CATE_MSE_eta = mean(CATE_SE_eta),
                CATE_MSE_e = mean(CATE_SE_e),
                CATE_MSE_e_param = mean(CATE_SE_e_param),
                CATE_MSE_kallus = mean(CATE_SE_kallus)),bias]

summ2[,':='(CATE_MSE_gsar1_ratio = CATE_MSE_gsar1/CATE_MSE_e,
            CATE_MSE_gsar2_ratio = CATE_MSE_gsar2/CATE_MSE_e,
            CATE_MSE_kappa1_ratio = CATE_MSE_kappa1/CATE_MSE_e,
            CATE_MSE_kappa2_ratio = CATE_MSE_kappa2/CATE_MSE_e,
            CATE_MSE_eta_ratio = CATE_MSE_eta/CATE_MSE_e,
            CATE_MSE_e_param_ratio = CATE_MSE_e_param/CATE_MSE_e,
            CATE_MSE_kallus_ratio = CATE_MSE_kallus/CATE_MSE_e)]


p2 <-  ggplot(summ2, aes(x = bias)) +
  geom_line(aes(y = CATE_MSE_eta_ratio,colour = "Ours"), size = 1.3) +
  # geom_line(aes(y =MSE_oberst_ratio ,colour = "Ours"), size = 1.3) +
  geom_line(aes(y = CATE_MSE_gsar1_ratio,colour =  "Green et al. 2005"),  size = 1.3) +
  # geom_line(aes(y = CATE_MSE_gsar2_ratio,colour = "gsar2"), size = 1.3) +
  geom_line(aes(y = CATE_MSE_kappa1_ratio,colour = "Rosenman et al. 2021"),  size = 1.3) +
  # geom_line(aes(y = CATE_MSE_kappa2_ratio,colour = "kappa2"),  size = 1.3) +
  geom_line(aes(y = CATE_MSE_e_param_ratio,colour = "Parametric method\nbenchmark"), linetype = "dashed", size = 1.3) +
  geom_line(aes(y = CATE_MSE_kallus_ratio,colour = "Kallus et al. 2018"),  size = 1.3) +
  scale_colour_manual("",
                      breaks = c("Ours", "Oberst", "Green et al. 2005", "gsar2","Rosenman et al. 2021","kappa2","Kallus et al. 2018","Parametric method\nbenchmark"),
                      values = c("#800000", "#002058", "#ed008c","#a79d96","#33b1e5","#0f7361","#f5ac3d","480080")) +
  xlab("bias") +
  geom_hline(yintercept = 1, linetype = "dashed",size = 1.3) +
  ggtitle("10 strata CATE overall MSE") +
  ylab("relative MSE") +
  theme_bw() +
  theme(text = element_text(family = "CMU Sans Serif"))
p2


ggsave(file.path(data_folder,"CATE_MSE_morevars.pdf"),
       plot = p2,
       device = cairo_pdf,
       dpi = 1280)


