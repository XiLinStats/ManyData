
# Start up  ---------------------------------------------------------------

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
library(doSNOW)
library(MASS)
library(labelled)
library(numDeriv)



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


# Data --------------------------------------------------------------------

data("STAR_Kallus")

unlabelled(STAR_Kallus)
STAR_Kallus <- as.data.table(STAR_Kallus)
STAR_Kallus[,`:=`(older = (birthyear < 1980 | (birthyear == 1980 & birthmonth <= 2)),
                  Sex = (gender == 2),
                  Race = (race == 2),
                  Trt = (g1classtype == 1))]
STAR_Kallus[, Y_perc := pnorm(Y,mean(STAR_Kallus$Y),sd(STAR_Kallus$Y))]

possUNC <- STAR_Kallus[race <= 2 & (g1surban == 1 | g1surban == 3)]


# Parameters  -------------------------------------------------------------

# Number of iterations (different sets of data)
n_boot <- 500

# List of eta
eta_list <- c(seq(0,0.2,0.04),seq(0.25,1,0.05))

# bias introduced by down-weighting the control individuals with low outcomes
# p = 0.3 means that we down weight the bottom 30% outcome
p_list <- c(0.0001,seq(0.025,0.5,0.025),0.6,0.8)

# sampling proportion from the experimental data (q = 0.1 means 10% sample)
q <- 0.1

# frugal parameterization
forms <- list(c(Sex ~ 1),
              list(Trt ~ Sex, older ~ 1),
              Y ~ Trt + older,
              ~ 1)
family <- list(c(5), c(5,5), 1, 1)

forms2 <- causl:::tidy_formulas(unlist(forms[-2]), kwd = "cop")
full_form <- causl:::merge_formulas(forms2)
msks <- causl:::masks(forms2,family = family[-2],full_form$wh)
vars <- causl:::lhs(unlist(forms2[-length(forms2)]))

# calculate ground truth
fit_gc <- fitCausal2(dat = possUNC, forms2, unlist(family[-2]))
ATE_gc <- fit_gc$pars$Y$beta[2]

# Simulation code ---------------------------------------------------------
run_sims <- function(){
  res <- matrix(0,1,5)

  for (p in p_list) {
    tryCatch({
      unlabelled(STAR_Kallus)
      STAR_Kallus[, Y_perc_cap := pmin(Y_perc,p)]
      wh <- sample(c(TRUE,FALSE), size = nrow(possUNC), replace = TRUE, prob = c(q, 1 - q))
      UNC <- possUNC[wh]

      CONF_t_all <- STAR_Kallus[race <= 2 & !(stdntid %in% UNC$stdntid) & (g1classtype == 1)]
      #Take a random sample based on capped percentile
      rand <- sample(x = 1:nrow(CONF_t_all),1000,replace = F,prob = CONF_t_all$Y_perc_cap)
      CONF_t <- CONF_t_all[rand]
      CONF_c <- STAR_Kallus[race <= 2 & !(stdntid %in% UNC$stdntid) & (g1classtype == 2)]
      CONF <- rbind(CONF_t,CONF_c)

      mm_exp <- model.matrix(full_form$formula, data = UNC)

      fit_exp <- fitCausal2(dat = UNC, forms2, unlist(family[-2]))
      fit_obs <- fitCausal2(dat = CONF, forms2, unlist(family[-2]))


      elpd <- -Inf
      opt_samples <- NA
      out_eta <- NA

      for (eta in eta_list) {
        tryCatch(
          {samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta, n_sample = 2000)
          elpd_sample <- calculate_elpd(samples,UNC[,vars,with = F],fit_exp$mm,msks,family = c(5,1), method = "WAIC",inCop = 1:2,
                                        family = unlist(family[-2][-length(family[-2])]),
                                        fam_cop = family[length(family)])

          if (elpd < elpd_sample$elpd_waic) {
            elpd <- elpd_sample$elpd_waic
            opt_samples <- samples
            opt_eta <- eta

          }
          }
          ,error = function(e) {cat("ERROR : p = ",p, "eta = ", eta, " skipped", "\n")})
      }
      if (sum(is.na(opt_samples)) >0 ){
        next
      }else{
        ATE_eta <- colMeans(opt_samples)[3]
        res <- rbind(res, c(p, opt_eta, fit_exp$par[3], fit_obs$par[3], ATE_eta))
      }
    },error = function(e) {cat("ERROR : p = ",p, " skipped", "\n")})
  }

  return(res[-1,])
}


# Run simulations in parallel ---------------------------------------------

start <- Sys.time()

print(paste0("start of parallelisation :",Sys.time()))

results1 <- foreach(
  i = 1:n_boot,
  .combine = 'rbind',
  .packages = c('mvtnorm','data.table','coda','ManyData','loo','causl','MASS','doParallel',"labelled","numDeriv")
) %dopar% {
  .GlobalEnv$family <- family
  .GlobalEnv$forms_obs <- forms
  .GlobalEnv$forms2 <- forms2
  .GlobalEnv$full_form <- full_form
  .GlobalEnv$msks <- msks
  .GlobalEnv$vars <- vars
  .GlobalEnv$calculate_elpd <- calculate_elpd
  .GlobalEnv$approx_posterior <- approx_posterior
  .GlobalEnv$fitCausal2 <- fitCausal2
  .GlobalEnv$STAR_Kallus <- STAR_Kallus
  .GlobalEnv$possUNC <- possUNC
  x <- run_sims()
  return(x)
}
print(Sys.time() - start)

results1 <- as.data.table(results1)
colnames(results1) <- c("p", "opt_eta","ATE_e", "ATE_o","ATE_eta")
saveRDS(results1,"data/STAR_res_v1.rds")

# Plotting ----------------------------------------------------------------

dat <- results1

summ <- dat[,.(count = .N,
               eta = mean(opt_eta),
               ATE_e = mean(ATE_e),
               ATE_o = mean(ATE_o),
               ATE_eta = mean(ATE_eta),
               opt_eta = mean(opt_eta),
               MSE_0 = mean((ATE_e - ATE_gc)^2),
               MSE_eta = mean((ATE_eta - ATE_gc)^2)),p][order(p)]

summ[,MSE_eta_ratio := MSE_eta/MSE_0]
summ[,bias := ATE_o - ATE_e]


p2 <- ggplot(dat[round(p,2) == 0.30]) +
  geom_density( aes(x = ATE_e, fill = "r"), alpha = 0.5) +
  geom_density( aes(x = ATE_eta, fill = "b"), alpha = 0.5) +
  geom_density( aes(x = ATE_o, fill = "c"), alpha = 0.5) +
  scale_fill_manual(name = "estimates", values = c( "r" = "#800000","b" = "#002058", "c" = "#0f7361"), labels = c("b" = "Combined data", "r" = "Unconfounded data only", "c" = "Confounded data only")) +
  ylab("density") +
  xlab("ATE estimates") +
  ggtitle("ATE distribution") +
  theme_bw() +
  theme(text = element_text( family = "CMU Sans Serif"))

ggsave(file.path(data_folder,"ATEdist.pdf"),
       plot = p2,
       device = cairo_pdf,
       dpi = 1280)
