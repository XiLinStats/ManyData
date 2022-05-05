# C = {sex}
# Z = {Health}
# X = {Treatment}
# Y = {Outcome}


# Observational study -----------------------------------------------------

# formula
forms_obs <- list(Z ~ C1,
                  list(X ~ Z*C1,C1~ 1),
                  Y ~ X*C1,
                  ~ X*C1 #copula
)


# parameters


pars_obs <- list(C1 = list(beta = 0),
                 Z = list(beta = c(0.2,0.6),phi = 1),
                 X = list(beta = c(0.5,0.6,0.1,0.4)),
                 Y = list(beta = c(0.6,0.8  ,0.2,0.3),phi = 1),
                 cop = list(beta = matrix(c(1,0,0,0),nrow = 4)))

set.seed(111)

data_obs <- causalSamp(2500,formulas = forms_obs,pars = pars_obs, family = list(1,c(5,5),1,1))

summary(glm(Z ~ C1, family = gaussian, data = data_obs))$coef
# Estimate Std. Error  t value     Pr(>|t|)
# (Intercept) 0.1805374 0.01996986  9.04049 2.188235e-19
# C1          0.5824543 0.02860437 20.36242 1.300400e-88

glmT <- glm(X ~ Z*C1, family = binomial, data = data_obs)
summary(glmT)

# Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.49217    0.04282  11.495  < 2e-16 ***
#   Z            0.63312    0.04645  13.630  < 2e-16 ***
#   C1           0.12968    0.06956   1.864 0.062264 .
# Z:C1         0.27249    0.07289   3.738 0.000185 ***

#IPW
ps <- fitted(glmT)
wt <- data_obs$X/ps + (1 - data_obs$X)/(1 - ps)

summary(svyglm(Y ~ X*C1, design=svydesign(~1, data=data_obs, weights= ~ wt)))$coef
# Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 0.6020416 0.03466124 17.369304 1.157523e-65
# X           0.7776120 0.04426501 17.567191 4.406576e-67
# C1          0.1194004 0.06815873  1.751799 7.986979e-02
# X:C1        0.3946419 0.07731291  5.104476 3.440156e-07


fit_obs<-fitCausal(dat = data_obs,
                   formulas = forms_obs[-2],
                   family = c(1,1,1),
                   sandwich = FALSE
)
# log-likelihood:  -13569.27
# Z ~ C1
# est.   s.e.
# (intercept) 0.181 0.0200
# C1          0.582 0.0286
# residual s.e.:  1.02 0.0205
#
# Y ~ X * C1
# est.   s.e.
# (intercept) 0.621 0.0314
# X           0.773 0.0376
# C1          0.114 0.0509
# X:C1        0.383 0.0576
# residual s.e.:  0.984 0.0201
#
# copula parameters:
#   cop ~ X + C1 + X:C1
# Error in dimnames(x) <- dn :
#   length of 'dimnames' [1] not equal to array extent

MLE_beta_obs<- fit_obs$pars$Y$beta
MLE_beta_obs[c(2,3)]<- fit_obs$pars$Y$beta[c(3,2)]


# Experimental study -----------------------------------------------------

# formula

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

set.seed(112)

data_exp <- causalSamp(250,formulas = forms_exp,pars = pars_exp, family = list(1,c(5,5),1,1))
mean(data_exp$X)
summary(glm(Z ~ C1, family = gaussian, data = data_exp))$coef

glmT <- glm(X ~ 1, family = binomial, data = data_exp)
summary(glmT)$coef
# Estimate Std. Error    z value  Pr(>|z|)
# (Intercept) -0.04800922  0.1265275 -0.3794369 0.7043635

#IPW
ps <- fitted(glmT)
wt <- data_exp$X/ps + (1 - data_exp$X)/(1 - ps)

summary(lm(Y ~ X*C1, data = data_exp, weights = wt))$coef
# Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 0.6745641  0.1347827 5.0048285 1.065787e-06
# X           0.9278107  0.1868262 4.9661705 1.277974e-06
# C1          0.2634456  0.1822594 1.4454432 1.496061e-01
# X:C1        0.2288059  0.2572040 0.8895892 3.745556e-01


fit_exp<-fitCausal(dat = data_exp,
          formulas = forms_exp[-2],
          family = c(1,1,1)
)
# log-likelihood:  -659.987
# Z ~ C1
# est.   s.e. sandwich
# (intercept) 0.425 0.0903   0.0924
# C1          0.425 0.1243   0.1240
# residual s.e.:  0.963 0.0871 0.0875
#
# Y ~ X * C1
# est.  s.e. sandwich
# (intercept) 0.663 0.122    0.122
# X           0.953 0.151    0.143
# C1          0.242 0.165    0.165
# X:C1        0.271 0.213    0.215
# residual s.e.:  1 0.0896 0.0773
#
# copula parameters:
#   cop ~ X + C1 + X:C1
# est.  s.e. sandwich
# (intercept)  1.1364 0.243    0.221
# C1           0.0818 0.320    0.296
# X            0.3605 0.325    0.319
# C1:X        -0.4911 0.447    0.424

MLE_exp<- fit_exp$pars$Y$beta
MLE_exp[c(2,3)]<- fit_exp$pars$Y$beta[c(3,2)]

# theta_sim_raw <- run_MH_MCMC(c(0.6,0.2,0.8,0.3), 500, data_obs, data_exp,
#                                         msks_obs,theta_obs,
#                                         msks_exp,theta_exp,MLE_beta_exp, 0.2*diag(4),0)
#


# MLE from obs + exp combined
data_comb <- rbind(data_obs,data_exp)
fit_comb<-fitCausal(dat = data_comb,
                   formulas = forms_exp[-2],
                   family = c(1,1,1),
                   sandwich = FALSE
)
# log-likelihood:  -14238.85
# Z ~ C1
# est.   s.e.
# (intercept) 0.191 0.0195
# C1          0.576 0.0279
# residual s.e.:  1.02 0.0199
#
# Y ~ X * C1
# est.   s.e.
# (intercept) 0.623 0.0304
# X           0.778 0.0364
# C1          0.133 0.0478
# X:C1        0.365 0.0546
# residual s.e.:  0.987 0.0196
#
# copula parameters:
#   cop ~ X + C1 + X:C1
# Error in dimnames(x) <- dn :
#   length of 'dimnames' [1] not equal to array extent
MLE_beta_comb<- fit_comb$pars$Y$beta
MLE_beta_comb[c(2,3)]<- fit_comb$pars$Y$beta[c(3,2)]
ATE_comb <- MLE_beta_comb[3] + MLE_beta_comb[4] * mean(data_comb$C1)


# Save --------------------------------------------------------------------

saveRDS(data_exp,file.path(data_dir,"/data_exp.rds"))
saveRDS(data_obs,file.path(data_dir,"/data_obs.rds"))


# Read --------------------------------------------------------------------
# data_exp <- readRDS(file.path(data_dir,"/data_exp.rds"))
# data_obs <- readRDS(file.path(data_dir,"/data_obs.rds"))
