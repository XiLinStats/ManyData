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
  q <- seq(0,1,1/K)[-1]
  table <- as.data.table(list(bin = q,
                              C1 = quantile(data_obs_b$C1,q)))
  table[bin == 1, C1 := Inf][,K := 1:.N]

  # stratify data
  setkey(table,C1)
  setkey(data_exp_b,C1)
  setkey(data_obs_b,C1)

  data_exp_b <- table[data_exp_b, roll = -Inf]
  data_obs_b <- table[data_obs_b, roll = -Inf]

  return(list(data_exp_b = data_exp_b, data_obs_b = data_obs_b))

}

# Params ------------------------------------------------------------------

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
family <- list(1, c(5,1), 1, 1)
forms_obs <- forms_exp <- list(Z ~ C1,
                               list(X ~ Z*C1 , C1 ~ 1),
                               Y ~ X*C1 ,
                               ~ X*C1)
forms_exp[[2]] <- list(X ~ 1, C1 ~ 1)

pars_exp <- list(C1 = list(beta = 0,phi = 1),
                 Z = list(beta = c(0.2, 0.6), phi = 1),
                 X = list(beta = c(0.5,0.6,0.1,0.4)),
                 Y = list(beta = c(0.6,0.8,0.2,0.3),phi = 1),
                 cop = list(beta = matrix(c(1,2.5,0,0),nrow = 4)))
pars_obs <- pars_exp
pars_exp$X <- list(beta = 0)

forms2 <- list(obs = causl:::tidy_formulas(forms_obs[-2], kwd = "cop"),
               exp = causl:::tidy_formulas(forms_exp[-2], kwd = "cop"))


full_form <- list(obs =  causl:::merge_formulas(forms2$obs),
                  exp =  causl:::merge_formulas(forms2$exp))

msks <- list(obs = ManyData:::masks(forms2$obs,family = family[-2],full_form$obs$wh),
             exp = ManyData:::masks(forms2$exp,family = family[-2],full_form$exp$wh))

theta <- list(obs = causl:::theta(pars = pars_obs, formulas = forms2$obs, full_form$obs, kwd = "cop"),
              exp = causl:::theta(pars = pars_exp, formulas = forms2$exp, full_form$exp, kwd = "cop"))
vars <- causl:::lhs(unlist(forms2$obs[1:2]))

# Simulation code ---------------------------------------------------------


# Number of iterations (different data)
# n_boot <- 44*5
n_boot <- 300

# List of eta
eta_list <- c(seq(0,0.2,0.04),seq(0.25,1,0.05))
# eta_list <- c(0.2)

bias_list <- c(0,0.5,0.75,1,1.25,1.5,1.75)

set.seed(111)

forms2 <- list(obs = causl:::tidy_formulas(forms_obs[-2], kwd = "cop"),
               exp = causl:::tidy_formulas(forms_exp[-2], kwd = "cop"))

full_form <- list(obs =  causl:::merge_formulas(forms2$obs),
                  exp =  causl:::merge_formulas(forms2$exp))

msks <- list(obs = ManyData:::masks(forms2$obs,family = family[-2],full_form$obs$wh),
             exp = ManyData:::masks(forms2$exp,family = family[-2],full_form$exp$wh))

theta <- list(obs = causl:::get_theta(pars = pars_obs, formulas = forms2$obs, full_form$obs, kwd = "cop"),
              exp = causl:::get_theta(pars = pars_exp, formulas = forms2$exp, full_form$exp, kwd = "cop"))

vars <- causl:::lhs(unlist(forms2$obs[1:2]))

func <- function(){


  result <- matrix(0,1,16)

  for (k in 1:length(bias_list)) {

    bias <- bias_list[k]

    coeff_UX <- bias
    coeff_UY <- bias

    family_d <- list(1, c(5,1), 1, 1)
    forms_obs_d <- forms_exp_d <- list(list(Z ~ C1),
                                       list(X ~ Z*C1 + U, C1 ~ 1),
                                       Y ~ X*C1 + U,
                                       ~ X*C1)
    forms_exp_d[[2]] <- list(X ~ 1, C1 ~ 1)

    pars_obs_d <- pars_exp_d <- list(C1 = list(beta = 0,phi = 1),
                                     Z = list(beta = c(0.2, 0.6), phi = 1),
                                     # Z2 = list(beta=c(0,0.4), phi=1-0.04^2),
                                     X = list(beta = c(0.5,0.6,0.1,coeff_UX,0.4)),
                                     Y = list(beta = c(0.6,0.8,0.2,coeff_UY,0.3),phi = 2),
                                     cop = list(beta = matrix(c(1,2.5,0,0),nrow = 4)))

    pars_exp_d$X <- list(beta = 0)

    data_obs_b <- data.frame(U = rbinom(n_o, size = 1, prob = 0.5))
    data_obs_b <- as.data.table(causl:::rfrugalParam(n = n_o, formulas = forms_obs_d, family = family_d, pars = pars_obs_d, dat = data_obs_b))

    data_exp_b <- data.frame(U = rbinom(n_e, size = 1, prob = 0.5))
    data_exp_b <- as.data.table(causl:::rfrugalParam(n = n_e, formulas = forms_exp_d, family = family_d, pars = pars_exp_d, dat = data_exp_b))


    # use the MLE fit in the experimental data to initialise the MCMC chain
    fit_exp <- fitCausal(dat = data_exp_b,
                         formulas = forms2$exp,
                         family = c(1,1,1)
    )

    fit_obs <- fitCausal(dat = data_obs_b, formulas = forms2$obs, family = c(1,1,1) )
    mm <- list(obs = fit_obs$mm, exp = fit_exp$mm)



    for (eta in eta_list){

      tryCatch(
        {samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta, n_sample = 2000)
        elpd_sample <- calculate_elpd(samples,data_exp_b[,vars,with = F],mm$exp, msks = msks$exp,method = "WAIC")
        }
        ,error = function(e) {cat("ERROR : Bias = ",bias, "eta = ", eta, " skipped", "\n")})

      if (is.na(samples)) {
        next
      }else{

        # Z grouping
        data_exp_b[,C1_group := 1*( C1 > 0)]
        C1_means <- c(mean(data_exp_b[C1_group == 1]$C1),mean(data_exp_b[C1_group == 0]$C1))


        # EDIT THIS
        result <- rbind(result,c(bias,eta,colMeans(samples[,c(1:7,9,11,12)]),mean(data_exp_b$C1),C1_means,elpd_sample$elpd_waic))
      }

    }
  }
  return(result[-1,])

}


# Run simulations in parallel ---------------------------------------------

start <- Sys.time()

print(paste0("start of parallelisation :",Sys.time()))

results1 <- foreach(
  i = 1:n_boot,
  .combine = 'rbind',
  .packages = c('mvtnorm','data.table','coda','ManyData','loo','causl','mvtnorm','MASS')
) %dopar% {
  .GlobalEnv$full_form <- full_form
  x <- func()
  return(x)
}
print(Sys.time() - start)
results1 <- as.data.table(results1)
colnames(results1) <- c("bias","eta","alpha_0","alpha_1","beta_0","beta_2","beta_1","beta_3","gamma_0","gamma_1","sigma_z","sigma_y",
                        "meanC1_exp","meanc1_1","meanc1_0", "ELPD_waic")
# print(results1)
saveRDS(results1,"data/results_v3_new.rds")



# Plotting ----------------------------------------------------------------

dat <- results1beta_1_true <- 0.8
beta_3_true <- 0.3


results1[,ind := ceiling((1:.N)/length(unique(result$eta)))]
eta_opt <- results1[ , .SD[which.max(ELPD_waic)], by = ind]
eta_opt[,mean(eta),bias]

create_plot <- function(result, title){

  result[,ATE := beta_1 + beta_3 * meanC1_exp]
  result[,CATE_C1 := beta_1 + beta_3 * meanc1_1]
  result[,CATE_C0 := beta_1 + beta_3 * meanc1_0]

  result[,ATE_true := beta_1_true + beta_3_true * meanC1_exp]
  result[,CATE_C1_true := beta_1_true + beta_3_true * meanc1_1]
  result[,CATE_C0_true := beta_1_true + beta_3_true * meanc1_0]


  summ <- result[,.(.N,
                    # eta = mean(eta),
                    ATE = mean(ATE),
                    elpd = mean(ELPD_waic),
                    MSE_C1 = mean((CATE_C1 - CATE_C1_true)^2),
                    MSE_C0 = mean((CATE_C0 - CATE_C0_true)^2),
                    MSE_ATE = mean((ATE - (beta_1_true + beta_3_true *meanC1_exp ))^2),
                    Bias = (mean(ATE - (beta_1_true + beta_3_true *meanC1_exp )))^2,
                    Variance = mean((ATE - (beta_1_true + beta_3_true *meanC1_exp )  - mean(ATE - (beta_1_true + beta_3_true *meanC1_exp )))^2),
                    Bias_C1= (mean(CATE_C1 - CATE_C1_true))^2,
                    Variance_C1 = mean((CATE_C1 - mean(CATE_C1))^2),
                    Bias_C0 = (mean(CATE_C0 - CATE_C0_true))^2,
                    Variance_C0 = mean((CATE_C0 - mean(CATE_C0))^2)

  ),

  ,.(eta,bias)]



  p1 <- ggplot(summ, aes(x = eta)) +
    geom_line(aes(y = -elpd), color = "#800000", size = 1.3) +
    ggtitle(title) +
    theme_bw() +
    theme(text = element_text(family = "CMU Sans Serif"))

  summ1 <- melt(summ, id.vars = c('eta'),
               measure.vars = c("Bias", "Variance"))

  p4 <- ggplot(summ1, aes(x = eta, y = value)) +
    geom_area(aes(fill = variable) , alpha = 0.9) +
    theme_bw() +
    labs(y = "MSE") +
    guides(fill = "none") +
    ggtitle("ATE") +
    theme(plot.title = element_text(size = 8)) +
    theme(text = element_text(family = "CMU Sans Serif")) +
    scale_fill_manual(values = c("#800000","#002058"))


  summ2 <- melt(summ, id.vars = c('eta'),
               measure.vars = c("Bias_C1", "Variance_C1"))

  p5 <- ggplot(summ2, aes(x = eta, y = value)) +
    geom_area(aes(fill = variable) , alpha = 0.9) +
    theme_bw() +
    labs(y = "MSE") +
    guides(fill = "none") +
    ggtitle("CATE for C > 0") +
    theme(plot.title = element_text(size = 8)) +
    theme(text = element_text(family = "CMU Sans Serif")) +
    scale_fill_manual(values = c("#800000","#002058"))

  return(list(p1,p4,p5,summ))

}

baseline_plots <- create_plot(result[ELPD_waic > -9e9 & bias == 0], title = "No influence")
baseline_plots[[1]] <- baseline_plots[[1]] + geom_vline(xintercept = 0.85, colour = "#f5ac3d", size = 1, linetype = "dotted")
baseline_plots[[2]] <- baseline_plots[[2]] + geom_vline(xintercept = 0.85, colour = "#f5ac3d", size = 1, linetype = "dotted")
baseline_plots[[3]] <- baseline_plots[[3]] + geom_vline(xintercept = 0.85, colour = "#f5ac3d", size = 1, linetype = "dotted")
small_plots <- create_plot(result[ELPD_waic > -9e9 & bias == 0.75], title = "Small influence")
medium_plots <- create_plot(result[ELPD_waic > -9e9 & bias == 1.25], title = "Large influence")

small_plots[[1]] <- small_plots[[1]]  + geom_vline(xintercept = 0.70, colour = "#f5ac3d", size = 1, linetype = "dotted")
medium_plots[[1]] <- medium_plots[[1]]  + geom_vline(xintercept = 0.18, colour = "#f5ac3d", size = 1, linetype = "dotted")

small_plots[[2]] <- small_plots[[2]]  + geom_vline(xintercept = 0.70, colour = "#f5ac3d", size = 1, linetype = "dotted")
medium_plots[[2]] <-  medium_plots[[2]]  + geom_vline(xintercept = 0.18, colour = "#f5ac3d", size = 1, linetype = "dotted")

small_plots[[3]] <- small_plots[[3]]  + geom_vline(xintercept = 0.70, colour = "#f5ac3d", size = 1, linetype = "dotted")
medium_plots[[3]] <-  medium_plots[[3]]  + geom_vline(xintercept = 0.18, colour = "#f5ac3d", size = 1, linetype = "dotted")


p <- grid.arrange(baseline_plots[[1]],baseline_plots[[2]], baseline_plots[[3]],
                  small_plots[[1]],small_plots[[2]],small_plots[[3]],
                  medium_plots[[1]],medium_plots[[2]],medium_plots[[3]],
                  ncol = 3)

ggsave(file.path(data_folder,"plots_new.pdf"),
       plot = p,
       device = cairo_pdf,
       dpi = 1280)
