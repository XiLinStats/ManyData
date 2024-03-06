##' Compute ELPD for particular value of learning rate
##'
##' @param eta learning rate
##' @param dat_e,dat_o data frames containing experimental and observational data
##' @param formulas formulas for distributions
##' @param family vector of family values
##' @param mode what should ELPD be evaluated upon (only \code{"exp"} works at the moment)
##' @param msks masks for observational and experimental data
##' @param full_formulas merged formulas for observational and experimental data
##' @param theta vector of parameters for observational and experimental data
##' @param start optional starting vector of parameters
##' @param mcmc_pars list of parameters for MCMC procedure
##'
##' @details The arguments \code{msks}, \code{full_formulas}, \code{theta} and
##' \code{start} are all optional, while \code{mcmc_pars} and \code{mode} have
##' defaults.
##'
##' @export
compute_elpd <- function (eta, dat_e, dat_o, formulas, family, mode="exp",
                          msks, full_formulas, theta, start, mcmc_pars=list()) {

  # get MCMC control parameters or use defaults
  mcp = list(n_iter = 1e3, n_burn=1e2, n_thin=5)
  matches = match(names(mcmc_pars), names(mcp))
  mcp[matches] = mcmc_pars[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in mcmc_pars not matched: ",
                                   paste(names(mcmc_pars[is.na(matches)]),
                                         sep = ", "))
  if (mcp$n_iter <= mcp$n_burn) stop("No samples will be retained")

  n_e <- nrow(dat_e)
  n_o <- nrow(dat_o)

  # use the MLE fit in the experimental data to initialise the MCMC chain
  if (missing(start)) {
    fit_exp <- fitCausal(dat_e, formulas[-2], family[-2])
    start <- fit_exp$par
  }

  if (missing(full_formulas)) {
    full_formulas <- list(obs=merge_formulas(formulas = formulas),
                          exp=merge_formulas(formulas = formulas))
  }
  if (missing(msks)) {
    msks <- list(obs=masks(formulas = formulas[-2],
                            family = family[-2],
                            wh = full_formulas$obs$wh[-2]),
                  exp=masks(formulas = formulas[-2],
                            family = family[-2],
                            wh = full_formulas$exp$wh[-2]))
  }
  if (missing(theta)) {
    theta <- list(obs=start,
                  exp=start)
    # theta <- list(obs=causl:::theta(pars = pars, formulas = formulas[-2],
    #                                 full_formulas$obs, kwd = "cop"),
    #               exp=causl:::theta(pars = pars, formulas = formulas[-2],
    #                                 full_formulas$exp, kwd = "cop"))
  }

  vars <- causl:::lhs(formulas[c(1,3)])

  theta_sim_raw <- run_MH_MCMC(start, mcmc_pars=mcp,
                               dat_o, dat_e,
                               msks$obs, theta$obs,
                               msks$exp, theta$exp,
                               start, 0.2*diag(length(start)),
                               eta,
                               full_formulas,
                               vars=vars)
  ni <- mcp$n_iter
  nb <- mcp$n_burn
  theta_sim_rtnd <- theta_sim_raw[seq_len(ni-nb) + nb, , drop=FALSE]
  if (mcp$n_thin > 1) theta_sim_rtnd <- theta_sim_rtnd[seq(to=nrow(theta_sim_rtnd), by=mcp$n_thin), , drop=FALSE]

  ## now estimate ELPD
  out <- matrix(NA, nrow(theta_sim_rtnd), n_e)

  np <- sum(msks$exp$beta_m > 0)
  mm_exp <- model.matrix(full_formulas$exp$formula, data = dat_e)

  # create the log likelihood matrix to calculate ELPD (Vehtari et al. 2017)
  for (i in seq_len(nrow(theta_sim_rtnd))) {

    # theta2 <- copy(theta_exp)
    theta2 <- theta_sim_rtnd[i,]
    msks2 <- msks$exp

    msks2$beta_m[msks2$beta_m > 0] <- theta2[seq_len(np)]
    msks2$phi_m[msks2$phi_m > 0] <- theta2[-seq_len(np)]
    out[i,] <- ManyData:::llC(dat_e[,vars],mm_exp,msks2$beta_m, phi = msks2$phi_m,inCop = c(1,2))
    # lst[[i]] <- as.vector(ll_i)
    # print(i)

  }

  ## Could add back in code to compute ESS etc here

  # mtrx <- do.call(rbind, lst)

  # calculate ELPD  (Vehtari et al. 2017)
  waic_eta <- loo::waic(out)
  # loo_eta <- loo::loo(out)

  waic_eta
}
