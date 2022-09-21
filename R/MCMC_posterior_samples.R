#' @export
masks <- function(formulas, family = rep(1, nc), wh, LHS)
{
  if (is.list(family)) {
    ncop <- length(family[length(family)])
    family <- unlist(family)
  }
  else {
    ncop <- 1
    # family <- family[-length(family)]
  }
  formulas <- unlist(formulas)
  nc <- length(formulas) + ncop - 1
  beta_m <- matrix(0, nrow = max(unlist(wh)), ncol = nc)
  phi_m <- numeric(length(family)-ncop)
  for (i in seq_along(phi_m)) {
    if (family[i] >= 1 && family[i] <= 3) {
      phi_m[i] <- 1
    }
    beta_m[wh[[i]], i] <- 1
  }

  cp <- length(phi_m) + 1
  for (i in seq_len(ncop)){
    beta_m[wh[[cp]], cp + i - 1] <- 1
  }
  return(list(beta_m = beta_m, phi_m = phi_m))
}



ll_sum <- function(masks, theta, mm, dat){

  # theta2 <- copy(theta)
  msks2 <- copy(masks)
  np <- sum(msks2$beta_m > 0)
  msks2$beta_m[masks$beta_m > 0] <- theta[seq_len(np)]
  msks2$phi_m[masks$phi_m > 0] <- theta[-seq_len(np)]
  ll_i <- causl:::ll(dat = dat[,c(1,5)], mm = mm, beta = msks2$beta_m, phi = msks2$phi_m, inCop = c(1,2),
                     fam_cop = 1, family = list(1,1), link = NULL, par2 = NULL,
                     useC = TRUE)
  return(sum(ll_i))
}

# proposalfunction <- function(param, sigma){
#
#   rmvnorm(1, param, sigma)
# }
#
# posterior <- function(simval, dat_obs, theta_obs,mm_obs, mask_obs,
#                       dat_exp, theta_exp,mm_exp,mask_exp,
#                       eta){
#
#   theta2_obs <- copy(theta_obs)
#   theta2_obs[6] <- simval[1]
#   theta2_obs[8] <- simval[2]
#
#   theta2_exp <- copy(theta_exp)
#   theta2_exp[6] <- simval[1]
#   theta2_exp[8] <- simval[2]
#
#   # likelihood from the observational data
#   ll_obs <- -causl:::nll2(
#     theta2_obs,
#     dat_obs[,c(1,5)],
#     mm_obs,
#     beta = mask_obs$beta_m,
#     phi = mask_obs$phi_m,
#     inCop = c(1,2),
#     fam_cop = 1,
#     family = list(1,1),
#     par2 = NULL,
#     useC = TRUE
#   )
#
#   # likelihood from the observational data
#   ll_exp <- -causl:::nll2(
#     theta2_exp,
#     dat_exp[,c(1,5)],
#     mm_exp,
#     beta = mask_exp$beta_m,
#     phi = mask_exp$phi_m,
#     inCop = c(1,2),
#     fam_cop = 1,
#     family = list(1,1),
#     par2 = NULL,
#     useC = TRUE
#   )
#
#   posterior <- dmvnorm(simval,mean = c(0,0),sigma = matrix(c(2,0,0,2), ncol = 2), log=TRUE) +
#     ll_exp +
#     eta * ll_obs
#
#   posterior
# }



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

  # Prepare input for the observational study
  ## get param masks
  # forms2_obs <- causl:::tidy_formulas(forms_obs[-2], kwd = "cop")
  # full_form_obs <- causl:::merge_formulas(forms2_obs)
  # wh_obs <- full_form_obs$wh
  # # LHS <- lhs(forms2[-length(forms2)])
  # msks_obs <- masks(forms_obs[-2],family = list(1,1,1),wh_obs)


  # theta_obs <- causl:::theta(pars = pars_obs, formulas = forms_obs[-2], full_form_obs, kwd = "cop")

  # Prepare input for the experimental study
  ## get param masks
  # forms2_exp <- causl:::tidy_formulas(forms_exp[-2], kwd = "cop")
  # full_form_exp <- causl:::merge_formulas(forms2_exp)
  # wh_exp <- full_form_exp$wh
  # # LHS <- lhs(forms2[-length(forms2)])
  # msks_exp <- masks(forms_exp[-2],family = list(1,1,1),wh_exp)


  # theta_exp <- causl:::theta(pars = pars_exp, formulas = forms_exp[-2], full_form_exp, kwd = "cop")
  if (missing(prop_cov)) {
    FI_exp <- ManyData:::ApproxFI_single(msks_exp, theta_exp, mm_exp, dat_exp, delta = 1e-4)
    FI_obs <- ManyData:::ApproxFI_single(msks_obs, theta_obs, mm_obs, dat_obs, delta = 1e-4)
    sigma <- 400* 1/(10*FI_exp+eta *FI_obs)

    # if (eta > 0.1){
    #
    # sigma <- 2 * solve(eta *FI_obs)
    #
    # } else{
    #   FI_exp <- ApproxFI(masks = msks_exp, theta = theta_exp, mm = mm_exp, dat = dat_exp, delta = 1e-4)
    #   sigma <- 2 * solve(FI_exp+eta *FI_obs)
    # }
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

  # ## fun MCMC
  # chain <- ManyData:::MCMCloop_C(n_iter = mcmc_pars$n_iter,
  #
  #                     init_val = startval,
  #                     # init_val = c(0,0,0,0),
  #                     sigma = sigma,
  #
  #                     inCop = seq_len(ncol(dat_exp)),
  #                     dat_exp = dat_exp,
  #                     theta_exp  = theta_exp,
  #                     mm_exp = mm_exp,
  #                     mask_exp = msks_exp,
  #
  #                     dat_obs = dat_obs,
  #                     theta_obs  = theta_obs,
  #                     mm_obs = mm_obs,
  #                     mask_obs = msks_obs,
  #
  #                     p_mu = p_mu,
  #                     p_sigma = p_sig,
  #                     eta = eta)

  return(chain)
}

MCMC_loop <- function (n_iter, n_burn = 0L, n_thin = 1L, init_val,
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

# run_MH_MCMC <- function(startval, iterations, dat_obs, dat_exp,
#                         msks_obs,theta_obs,
#                         msks_exp,theta_exp,
#                         p_mu, p_sig,eta){
#
#   # Prepare input for the observational study
#   ## get param masks
#   # forms2_obs <- causl:::tidy_formulas(forms_obs[-2], kwd = "cop")
#   # full_form_obs <- causl:::merge_formulas(forms2_obs)
#   # wh_obs <- full_form_obs$wh
#   # # LHS <- lhs(forms2[-length(forms2)])
#   # msks_obs <- masks(forms_obs[-2],family = list(1,1,1),wh_obs)
#
#   ## get model matrix
#   mm_obs <- model.matrix(full_form_obs$formula, data = dat_obs)
#
#   # theta_obs <- causl:::theta(pars = pars_obs, formulas = forms_obs[-2], full_form_obs, kwd = "cop")
#
#   # Prepare input for the experimental study
#   ## get param masks
#   # forms2_exp <- causl:::tidy_formulas(forms_exp[-2], kwd = "cop")
#   # full_form_exp <- causl:::merge_formulas(forms2_exp)
#   # wh_exp <- full_form_exp$wh
#   # # LHS <- lhs(forms2[-length(forms2)])
#   # msks_exp <- masks(forms_exp[-2],family = list(1,1,1),wh_exp)
#
#   ## get model matrix
#   mm_exp <- model.matrix(full_form_exp$formula, data = dat_exp)
#
#   # theta_exp <- causl:::theta(pars = pars_exp, formulas = forms_exp[-2], full_form_exp, kwd = "cop")
#
#
#   # covariance matrix for proposal function
#   FI_exp <- ManyData:::ApproxFI_all(masks = msks_exp, theta = theta_exp, mm = mm_exp, dat = dat_exp, delta = 1e-4)
#   FI_obs <- ManyData:::ApproxFI_all(masks = msks_obs, theta = theta_obs, mm = mm_obs, dat = dat_obs, delta = 1e-4)
#
#   # if (eta > 0.1){
#   #
#   # sigma <- 2 * solve(eta *FI_obs)
#   #
#   # } else{
#   #   FI_exp <- ApproxFI(masks = msks_exp, theta = theta_exp, mm = mm_exp, dat = dat_exp, delta = 1e-4)
#   #   sigma <- 2 * solve(FI_exp+eta *FI_obs)
#   # }
#   sigma <- 2 * solve(FI_exp + eta *FI_obs)
#
#   # chain <- matrix(rep(0,iterations*2), ncol = 2)
#   # chain[1,] <- startval
#   #
#   # for (i in 2:iterations) {
#   #
#   #   currentval = chain[i - 1,]
#   #   Y = proposalfunction(currentval, sigma)
#   #
#   #   alpha = exp(posterior(Y, dat_obs = dat_obs, theta_obs = theta_obs, mm_obs = mm_obs, mask_obs = msks_obs,
#   #                         dat_exp = dat_exp, theta_exp = theta_exp,mm_exp = mm_exp ,mask_exp = msks_exp,
#   #                         eta)
#   #               - posterior(currentval, dat_obs = dat_obs,theta_obs = theta_obs,  mm_obs = mm_obs, mask_obs = msks_obs,
#   #                           dat_exp = dat_exp, theta_exp = theta_exp,mm_exp = mm_exp ,mask_exp = msks_exp,
#   #                           eta))
#   #
#   #   if (runif(1) < alpha) {
#   #     chain[i,] = Y
#   #   }else{
#   #     chain[i,] = currentval
#   #   }
#   #
#   #   # printCount(i)
#   # }
#
#   chain <- ManyData:::MCMCloop_C(n_iter = iterations,
#
#                                  init_val = startval,
#                                  # init_val = c(0,0,0,0),
#                                  sigma = sigma,
#
#                                  inCop = c(1,2),
#                                  dat_exp = dat_exp[,c("Z","Y")],
#                                  theta_exp  = theta_exp,
#                                  mm_exp = mm_exp,
#                                  mask_exp = msks_exp,
#
#                                  dat_obs = dat_obs[,c("Z","Y")],
#                                  theta_obs  = theta_obs,
#                                  mm_obs = mm_obs,
#                                  mask_obs = msks_obs,
#
#                                  p_mu = p_mu,
#                                  p_sigma = p_sig,
#                                  eta = eta)
#
#   return(chain)
#
# }
