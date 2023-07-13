#' Simulate posterior samples using MCMC
#'
#' @param startval Initial values.
#' @param n_iter Number of MCMC iterations.
#' @param dat_obs Observational dataset.
#' @param dat_exp Experimental dataset.
#' @param msks_obs,msks_exp masks for parameter vectors with observational and experimental distributions
#' @param theta_obs,theta_exp masks for parameter vectors with observational and experimental distributions
#' @param p_mu,p_sigma vector and matrix for Gaussian proposal distribution
#' @param eta Learning rate on the observational data.
# @param full_formulas list containing full formulae for \code{exp} and \code{obs} data
#' @param prop_cov optional covariance for proposal
#'
#' @details Proposal is currently an arbitrary multiple of univariate second
#' derivatives (10).
#'
#' @return A table of posterior samples.
#' @examples theta_sim_eta05 <- run_MH_MCMC(c(1,0.5), 1000, dat_obs = dat_obs, dat_exp = dat_exp, eta = 0.5)
#' @export
#'
#'
run_MH_MCMC <- function(startval, mcmc_pars,
                        dat_obs, dat_exp,
                        msks_obs, theta_obs, mm_obs,
                        msks_exp, theta_exp, mm_exp,
                        p_mu=rep(0,length(startval)),
                        p_sig=diag(length(startval)),
                        family, link,
                        eta, prop_cov) {
  if (missing(prop_cov)) {
    FI_exp <- ManyData:::ApproxFI_single(msks_exp, theta_exp, mm_exp, dat_exp, delta = 1e-4)
    FI_obs <- ManyData:::ApproxFI_single(msks_obs, theta_obs, mm_obs, dat_obs, delta = 1e-4)
    sigma <- 400* 1/(10*FI_exp+eta *FI_obs)
  }
  else sigma <- solve(prop_cov)

  chain <- MCMC_loop(n_iter = mcmc_pars$n_iter, n_burn = mcmc_pars$n_burn,
                     n_thin = mcmc_pars$n_thin, init_val = startval,
                     dat_e = dat_exp, theta_e = theta_exp, mm_e = mm_exp, mask_e = msks_exp,
                     dat_o = dat_obs, theta_o = theta_obs, mm_o = mm_obs, mask_o = msks_obs,
                     family = family,
                     eta = eta, prop_sigma = sigma)

  acc_rt <- sum(unique(chain))/nrow(chain)
  attr(chain, "AR") <- acc_rt

  return(chain)
}

#' @export
MCMC_loop <- function(n_iter, n_burn = 0L, n_thin = 1L, init_val,
                       dat_e, theta_e, mm_e, mask_e,
                       dat_o, theta_o, mm_o, mask_o,
                       family, link,
                       eta, prop_sigma) {
  if (missing(prop_sigma)) {
    stop("Must provide covariance matrix for proposal")
  }

  theta <- init_val

  curr_ll <- - causl:::nll2(theta, dat_e, mm_e, mask_e$beta_m, mask_e$phi_m, seq_along(mask_e$phi_m),
                  fam_cop=1, family, link, useC=TRUE) +
    - eta*causl:::nll2(theta, dat_o, mm_o, mask_o$beta_m, mask_o$phi_m, seq_along(mask_o$phi_m),
             fam_cop=1, family, link, useC=TRUE)

  out <- matrix(NA, nrow=ceiling((n_iter-n_burn)/n_thin), ncol=length(theta))

  rec <- 0L

  for (i in seq_len(n_iter)) {
    rje::printCount(i)
    mv <- mvtnorm::rmvnorm(n=1, sigma=prop_sigma/100)[1,]
    theta_prop <- theta + mv

    prop_ll <- - causl:::nll2(theta_prop, dat_e, mm_e, mask_e$beta_m, mask_e$phi_m, seq_along(mask_e$phi_m),
                    fam_cop=1, family, link, useC=TRUE) +
      - eta*causl:::nll2(theta_prop, dat_o, mm_o, mask_o$beta_m, mask_o$phi_m, seq_along(mask_o$phi_m),
               fam_cop=1, family, link, useC=TRUE)

    ## generate a log-uniform r.v.
    lU <- -rexp(1)

    ## check if new value should be accepted
    if (lU < prop_ll - curr_ll) {
      # print(theta)
      theta <- theta_prop
      curr_ll <- prop_ll
    }

    ## only record selected iterations
    if (i > n_burn && ((i - n_burn - 1) %% n_thin == 0)) {
      rec <- rec + 1
      out[rec,] <- theta
    }
  }

  out
}


#' @export
MCMC_loop_Gibbs <- function(n_iter, n_burn = 0L, n_thin = 1L, init_val,
                            dat_e, theta_e, mm_e, mask_e,
                            dat_o, theta_o, mm_o, mask_o,
                            family, link,
                            eta, prop_sigma) {
  if (missing(prop_sigma)) {
    stop("Must provide covariance matrix for proposal")
  }

  theta <- init_val

  curr_ll <- -causl:::nll2(theta, dat_e, mm_e, mask_e$beta_m, mask_e$phi_m, seq_along(mask_e$phi_m),
                           fam_cop = 1, family, link, useC = TRUE) +
    -causl:::nll2(theta, dat_o, mm_o, mask_o$beta_m, mask_o$phi_m, seq_along(mask_o$phi_m),
                  fam_cop = 1, family, link, useC = TRUE)

  out <- matrix(NA, nrow = ceiling((n_iter - n_burn)/n_thin), ncol = length(theta))

  rec <- 0L

  for (i in seq_len(n_iter)) {
    # rje::printCount(i)
    for (j in 1:length(init_val)) {

      theta_prop <- theta

      theta_prop[j] <- rnorm(1,theta[j],sqrt(prop_sigma[j,j]))

      prop_ll <- -causl:::nll2(theta_prop, dat_e, mm_e, mask_e$beta_m, mask_e$phi_m, seq_along(mask_e$phi_m),
                                fam_cop = 1, family, link, useC = TRUE) +
        -eta*causl:::nll2(theta_prop, dat_o, mm_o, mask_o$beta_m, mask_o$phi_m, seq_along(mask_o$phi_m),
                          fam_cop = 1, family, link, useC = TRUE)

      ## generate a log-uniform r.v.
      lU <- runif(1)

      ## check if new value should be accepted
      if (lU < exp(prop_ll - curr_ll)) {
        # print(theta)
        theta <- theta_prop
        curr_ll <- prop_ll
      }

    }

    ## only record selected iterations
    if (i > n_burn && ((i - n_burn - 1) %% n_thin == 0)) {
      rec <- rec + 1
      out[rec,] <- c(theta)
    }


  }

  out
}

#' Simulate posterior samples using Gibbs-MH MCMC
#'
#' @param startval Initial values.
#' @param n_iter Number of MCMC iterations.
#' @param dat_obs Observational dataset.
#' @param dat_exp Experimental dataset.
#' @param msks_obs,msks_exp masks for parameter vectors with observational and experimental distributions
#' @param theta_obs,theta_exp masks for parameter vectors with observational and experimental distributions
#' @param p_mu,p_sigma vector and matrix for Gaussian proposal distribution
#' @param eta Learning rate on the observational data.
# @param full_formulas list containing full formulae for \code{exp} and \code{obs} data
#' @param prop_cov optional covariance for proposal
#'
#' @details Proposal is currently an arbitrary multiple of univariate second
#' derivatives (10).
#'
#' @return A table of posterior samples.
#' @examples theta_sim_eta05 <- run_MH_MCMC(c(1,0.5), 1000, dat_obs = dat_obs, dat_exp = dat_exp, eta = 0.5)
#' @export
#'
#'
run_MH_MCMC_Gibbs <- function(startval, mcmc_pars,
                              dat_obs, dat_exp,
                              msks_obs, theta_obs, mm_obs,
                              msks_exp, theta_exp, mm_exp,
                              p_mu=rep(0,length(startval)),
                              p_sig=diag(length(startval)),
                              family, link,
                              eta, prop_cov) {

  if (missing(prop_cov)) {
    FI_exp <- ManyData:::ApproxFI_all(msks_exp, theta_exp, mm_exp, dat_exp, delta = 1e-4)
    FI_obs <- ManyData:::ApproxFI_all(msks_obs, theta_obs, mm_obs, dat_obs, delta = 1e-4)
    sigma <- 400* 1/(10*FI_exp+eta *FI_obs)

  }
  else sigma <- solve(prop_cov)

  chain <- MCMC_loop_Gibbs(n_iter = mcmc_pars$n_iter, n_burn = mcmc_pars$n_burn,
                           n_thin = mcmc_pars$n_thin, init_val = startval,
                           dat_e = dat_exp, theta_e = theta_exp, mm_e = mm_exp, mask_e = msks_exp,
                           dat_o = dat_obs, theta_o = theta_obs, mm_o = mm_obs, mask_o = msks_obs,
                           family = family,
                           eta = eta, prop_sigma = sigma)

  # acc_rt <- sum(unique(chain))/nrow(chain)
  # attr(chain, "AR") <- acc_rt

  return(chain)
}


#' Calculate the estimated ELPD
#'
#' @param samples posterior samples
#' @param data data to evaluate ELPD on, recommend to be the randomized data
#' @param mm model matrix of the data to evaluate ELPD on
#' @param msks masks for data, containing a beta matrix and a phi vector
#' @param method the method to estimate ELPD which takes values of "WAIC" or "LOO". See details.
#' @param inCop columns to include in the copula
#' @param family families of variable distributions
#' @param fam_cop copula family
#' @return Results from loo::waic() or loo::loo()
#' @examples elpd_sample <- calculate_elpd(samples,data_exp_b,mm$exp, msks = msks$exp,method = "WAIC")
#' @export


calculate_elpd <- function(samples,data,mm,msks,method = "WAIC",inCop, family, fam_cop){
  if (!(method %in% c("WAIC","LOO"))) {
    stop("Method has to be either WAIC LOO.")
  }
  
  lst <- list()
  
  for (i in 1:nrow(samples)) {
    
    # theta2 <- copy(theta_exp)
    theta2 <- samples[i,]
    msks2 <- copy(msks)
    np <- sum(msks$beta_m > 0)
    msks2$beta_m[msks$beta_m > 0] <- theta2[seq_len(np)]
    msks2$phi_m[msks$phi_m > 0] <- theta2[-seq_len(np)]
    
    ll_i <- causl:::ll(data, mm, msks2$beta_m, phi = msks2$phi_m, 
                       inCop = inCop,family = family, fam_cop = fam_cop)
    lst[[i]] <- as.vector(ll_i)
  }
  
  mtrx <- do.call(rbind, lst)
  
  if (method == "WAIC") {
    elpd <- loo::waic(mtrx)
  }
  if (method == "LOO") {
    elpd <- loo::loo(mtrx)
  }
  return(elpd)
  
}


#' Simulate posterior samples through normal approximation
#'
#' @param fit_exp The fitCausal output from the experimental data.
#' @param fit_obs The fitCausal output from the observational data.
#' @param eta learning rate eta.
#' @param n_sample Number of posterior samples to be generated from the normal distribution.
#' @return A matrix of n_sample by length(theta) matrix with one sample in each row.
#' @examples samples <- approx_posterior(fit_exp = fit_exp,  fit_obs, eta = eta, n_sample = 5000)
#' @export


approx_posterior <- function(fit_exp, fit_obs, eta, n_sample){

  if (eta < 0) {
    stop("Must provide a positive eta")
  }

  if (eta > 1) {
    warning("We recommend setting eta between 0 and 1")
  }

  V <- n_e * fit_exp$sandwich

  if (!any(is.na(V)) && rcond(V) > 1e-16) {
    invV <- solve.default(V)
  }
  else {
    invV <- tryCatch(MASS::ginv(V), error = function(e) {cat("ERROR : V is singular", "\n")})
  }

  W <- n_o * fit_obs$sandwich

  if (!any(is.na(W)) && rcond(W) > 1e-16) {
    invW <- solve.default(W)
  }
  else {
    invW <- tryCatch(MASS::ginv(W), error = function(e) {cat("ERROR : W is singular", "\n")})
  }

  theta_0 <- fit_exp$par
  theta_inf <- fit_obs$par

  VW <- n_e * invV + eta * n_o * invW

  if (!any(is.na(VW)) && rcond(VW) > 1e-16) {
    invvW <- solve.default(VW)
  }
  else {
    invvW <- tryCatch(MASS::ginv(VW), error = function(e) {cat("ERROR : VW is singular", "\n")})
  }


  theta_hat <- invvW %*% (n_e * invV %*% theta_0 + eta * n_o * invW %*% theta_inf)

  FI <- fit_exp$FI + eta * fit_obs$FI

  if (!any(is.na(FI)) && rcond(FI) > 1e-16) {
    invFI <- solve.default(FI)
  }
  else {
    invFI <- tryCatch(MASS::ginv(FI), error = function(e) {cat("ERROR : FI is singular", "\n")})
  }

  samples <- mvrnorm(n_sample,theta_hat,invFI)

  return(samples)
}
