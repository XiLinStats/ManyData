
# Input -------------------------------------------------------------------

dat  <- data_exp
formulas <- forms_exp[-2]
family <- c(1,1,1)
par2 <- NULL
pars<-pars_exp


forms <- causl:::tidy_formulas(formulas, kwd = "cop")
fam_cop <- last(family)
LHS <- causl:::lhs(forms[-length(forms)])
full_form <- causl:::merge_formulas(forms)
wh <- full_form$wh
mm <- model.matrix(full_form$formula, data = dat)
theta <- causl:::theta(pars = pars, formulas = formulas, full_form, kwd = "cop")
fix_par <- theta



# Hack nll2 ---------------------------------------------------------------

nll3 <- function (par, fix_par, dat, mm, beta, phi, inCop, fam_cop = 1, family = rep(1,
                                                                      nc), link, par2 = NULL, useC = TRUE)
{

  theta <- copy(fix_par)
  theta[3] <- par[1]
  theta[4] <- par[2]
  theta[5] <- par[3]
  theta[6] <- par[4]
  np <- sum(beta > 0)
  beta[beta > 0] <- theta[seq_len(np)]
  phi[phi > 0] <- theta[-seq_len(np)]
  -sum(causl:::ll(dat, mm = mm, beta = beta, phi = phi, inCop = inCop,
          fam_cop = fam_cop, family = family, link = link, par2 = par2,
          useC = useC))
}



# FitCausal ---------------------------------------------------------------
con <- list(method = "BFGS", newton = FALSE, cop = "cop",
            trace = 0, fnscale = 1, maxit = 10000L, abstol = -Inf,
            reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
            gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
            lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10)

# length(theta_exp)
# length(theta_st)
#
beta_start2 <- causl:::initializeParams2(dat, formulas = forms,
                                 family = family, full_form = full_form, kwd = "cop")
#
#
# theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0],
#               beta_start2$phi[beta_start2$phi_m > 0])
other_args2 <- list(fix_par = fix_par,dat = dat[, LHS, drop = FALSE], mm = mm,
                    beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                    inCop = seq_along(LHS), fam_cop = fam_cop, fam = family[-length(family)],
                    par2 = par2, useC = T)
# do.call(causl:::nll2, c(list(theta = theta_exp), other_args2))
# # [1] 687.5311
# do.call(nll3, c(list(par = c(0.2,0.8)), other_args2))
# [1] 687.5311
# out$par <- c(0.4,0.4,0.4,0.4)
init <- c(0.4,0.4,0.4,0.4)
con$maxit <- 5000
out <- do.call(optim, c(list(fn = nll3, par = init),
                        other_args2, list(method = "Nelder-Mead", control = con)))

MLE_fix_exp <- out$par

# maxit <- con$maxit
# conv <- FALSE
# while (!conv) {
#   con$maxit <- 5000
#   out <- do.call(optim, c(list(fn = nll3, par = out$par),
#                           other_args2, list(method = "Nelder-Mead", control = con)))
#   con$maxit <- max(maxit - 5000, 100)
#   out2 <- tryCatch(do.call(optim, c(list(fn = nll3, par = out$par),
#                                     other_args2, list(method = "BFGS", control = con))),
#                    warning = function(e) NA, error = function(e) NA)
#   if (!isTRUE(is.na(out2))) {
#     out <- out2
#     conv <- TRUE
#   }
#   else out2 <- list(par = out$par)
# }
# curr_val = out$value
# out <- do.call(optim, c(list(fn = nll3, par = out$par),
#                         other_args2, list(method = "Nelder-Mead", control = con)))

# out2 <- tryCatch(do.call(optim, c(list(fn = nll3, par = out$par),
#                                   other_args2, list(method = "BFGS", control = con))),
#                  warning = function(e) NA, error = function(e) NA)
#
#

# Input -------------------------------------------------------------------

dat  <- data_obs
formulas <- forms_exp[-2]
family <- c(1,1,1)
par2 <- NULL
pars<-forms_exp


forms <- causl:::tidy_formulas(formulas, kwd = "cop")
fam_cop <- last(family)
LHS <- causl:::lhs(forms[-length(forms)])
full_form <- causl:::merge_formulas(forms)
wh <- full_form$wh
mm <- model.matrix(full_form$formula, data = dat)
theta <- causl:::theta(pars = pars, formulas = formulas, full_form, kwd = "cop")
fix_par <- theta




# FitCausal ---------------------------------------------------------------
con <- list(method = "BFGS", newton = FALSE, cop = "cop",
            trace = 0, fnscale = 1, maxit = 10000L, abstol = -Inf,
            reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
            gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
            lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10)

# length(theta_exp)
# length(theta_st)
#
beta_start2 <- causl:::initializeParams2(dat, formulas = forms,
                                         family = family, full_form = full_form, kwd = "cop")
#
#
# theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0],
#               beta_start2$phi[beta_start2$phi_m > 0])
other_args2 <- list(fix_par = fix_par,dat = dat[, LHS, drop = FALSE], mm = mm,
                    beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                    inCop = seq_along(LHS), fam_cop = fam_cop, fam = family[-length(family)],
                    par2 = par2, useC = T)
# do.call(causl:::nll2, c(list(theta = theta_exp), other_args2))
# # [1] 687.5311
# do.call(nll3, c(list(par = c(0.2,0.8)), other_args2))
# [1] 687.5311
# out$par <- c(0.4,0.4,0.4,0.4)
init <- c(0.4,0.4,0.4,0.4)
con$maxit <- 5000
out <- do.call(optim, c(list(fn = nll3, par = init),
                        other_args2, list(method = "Nelder-Mead", control = con)))
MLE_fix_obs <- out$par


# Input -------------------------------------------------------------------

dat  <- data_comb
formulas <- forms_exp[-2]
family <- c(1,1,1)
par2 <- NULL
pars<-forms_exp


forms <- causl:::tidy_formulas(formulas, kwd = "cop")
fam_cop <- last(family)
LHS <- causl:::lhs(forms[-length(forms)])
full_form <- causl:::merge_formulas(forms)
wh <- full_form$wh
mm <- model.matrix(full_form$formula, data = dat)
theta <- causl:::theta(pars = pars, formulas = formulas, full_form, kwd = "cop")
fix_par <- theta




# FitCausal ---------------------------------------------------------------
con <- list(method = "BFGS", newton = FALSE, cop = "cop",
            trace = 0, fnscale = 1, maxit = 10000L, abstol = -Inf,
            reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
            gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
            lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10)

# length(theta_exp)
# length(theta_st)
#
beta_start2 <- causl:::initializeParams2(dat, formulas = forms,
                                         family = family, full_form = full_form, kwd = "cop")
#
#
# theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0],
#               beta_start2$phi[beta_start2$phi_m > 0])
other_args2 <- list(fix_par = fix_par,dat = dat[, LHS, drop = FALSE], mm = mm,
                    beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                    inCop = seq_along(LHS), fam_cop = fam_cop, fam = family[-length(family)],
                    par2 = par2, useC = T)
# do.call(causl:::nll2, c(list(theta = theta_exp), other_args2))
# # [1] 687.5311
# do.call(nll3, c(list(par = c(0.2,0.8)), other_args2))
# [1] 687.5311
# out$par <- c(0.4,0.4,0.4,0.4)
init <- c(0.4,0.4,0.4,0.4)
con$maxit <- 5000
out <- do.call(optim, c(list(fn = nll3, par = init),
                        other_args2, list(method = "Nelder-Mead", control = con)))
MLE_fix_comb <- out$par

