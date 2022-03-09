# C = {sex,age}
# Z = {Health}
# X = {Treatment}
# Y = {Outcome}


# Observational study -----------------------------------------------------

# formula

forms_obs <- list(Z ~ C1*C2,
              list(X ~ C1*C2*Z,C1~ 1, C2 ~ 1),
              Y ~ X*C1*C2,
              ~ X*C1*C2  #copula
              )


# parameters
pars_obs <- list(C1 = list(beta = 0.5),
             C2 = list(beta = 2,phi = 0.5),
             Z = list(beta = c(0.2,0.6,0.9,0.3),phi = 1),
             X = list(beta = c(0.5,0.4,0.2,0.3,0,0,0,0)),
             Y = list(beta = c(0.6,0.8 ,0.2,0.3,0,0,0,0),phi = 1),
             cop = list(beta = matrix(c(1,0,0,0,0,0,0,0),nrow = 8)))

set.seed(112)

data_obs <- causalSamp(1e4,formulas = forms_obs,pars = pars_obs, family = list(1,c(5,5,1),1,1))


summary(glm(Z ~ C1*C2, family = gaussian, data = data_obs))$coef
# Estimate Std. Error   t value      Pr(>|t|)
# (Intercept) 0.2259576 0.04848255  4.660597  3.193711e-06
# C1          0.5780516 0.06186441  9.343848  1.127684e-20
# C2          0.8982749 0.02276577 39.457259 1.391558e-316
# C1:C2       0.2957819 0.02905811 10.178982  3.229883e-24

mean(data_obs$X)
# [1] 0.8605
glmT <- glm(X ~ C1*C2*Z, family = binomial, data = data_obs)
summary(glmT)$coef
# Estimate Std. Error     z value     Pr(>|z|)
# (Intercept)  0.250875353 0.17131285  1.46442807 1.430770e-01
# C1           0.677893720 0.26950185  2.51535833 1.189115e-02
# C2           0.315135880 0.10388108  3.03362143 2.416374e-03
# Z            0.365787322 0.09270315  3.94579162 7.953677e-05
# C1:C2       -0.145840059 0.16779583 -0.86915189 3.847641e-01
# C1:Z         0.012460359 0.12453320  0.10005652 9.202995e-01
# C2:Z        -0.036687310 0.04440212 -0.82625122 4.086616e-01
# C1:C2:Z      0.001835339 0.05907600  0.03106742 9.752158e-01

#IPW
ps <- fitted(glmT)
wt <- data_obs$X/ps + (1 - data_obs$X)/(1 - ps)
summary(lm(Y ~ X+C1+C2, data = data_obs, weights = wt))$coef

# Estimate Std. Error  t value      Pr(>|t|)
# (Intercept) 0.5897961 0.03378515 17.45726  2.998519e-67
# X           0.7768122 0.01979112 39.25054 1.601510e-313
# C1          0.2114160 0.02040021 10.36342  4.869363e-25
# C2          0.3038287 0.01397440 21.74180 1.905894e-102

summary(svyglm(Y ~ X+C1+C2, design=svydesign(~1, data=data_obs, weights= ~ wt)))$coef
# Estimate Std. Error   t value      Pr(>|t|)
# (Intercept) 0.5897961 0.05008788 11.775227  8.492282e-32
# X           0.7768122 0.03149059 24.668078 1.780820e-130
# C1          0.2114160 0.03032327  6.972069  3.320591e-12
# C2          0.3038287 0.02259345 13.447646  7.195448e-41

# fitCausal(dat = data_obs,
#           formulas = forms_obs[-2],
#           family = c(1,1,1)
# )


# Experimental study -----------------------------------------------------

# formula

forms_exp <- list(Z ~ C1*C2,
              list(X ~ C1*C2, C1~ 1, C2 ~ 1),
              Y ~ X*C1*C2,
              ~ X*C1*C2  #copula
)


# parameters
pars_exp <- list(C1 = list(beta = 0.5),
             C2 = list(beta = 2,phi = 1),
             Z = list(beta = c(0.2,0.6,0.9,0.3),phi = 1),
             X = list(beta = c(0.3,0.8,0.2,0.1)),
             Y = list(beta = c(0.6,0.8 ,0.2,0.3,0,0,0,0),phi = 1),
             cop = list(beta = matrix(c(1,0,0,0,0,0,0,0),nrow = 8)))

data_exp <- setDT(causalSamp(250,formulas = forms_exp,pars = pars_exp, family = list(1,c(5,5,1),1,1)))

summary(glm(Z ~ C1*C2, family = gaussian, data = data_exp))$coef
# Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 0.1171791 0.23329008 0.5022894 6.159133e-01
# C1          0.7990946 0.29409672 2.7171148 7.053220e-03
# C2          0.9429862 0.09602417 9.8203005 2.033433e-19
# C1:C2       0.2008297 0.12649530 1.5876457 1.136509e-01
# Say only those with health score above 0.6 are slected into the RCT

# data_exp2 <- data_exp[Z > 0.6]
# summary(data_exp2$Z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.603   1.024   1.510   1.614   2.046   4.455

mean(data_exp$X)
# [1] 0.706
glmT <- glm(X ~ C1*C2, family = binomial, data = data_exp)
summary(glmT)$coef
# Estimate Std. Error     z value  Pr(>|z|)
# (Intercept)  0.64285701  0.5112230  1.25748833 0.2085769
# C1           0.52295206  0.6758714  0.77374494 0.4390816
# C2           0.21136305  0.2213110  0.95504984 0.3395525
# C1:C2       -0.01859143  0.3073900 -0.06048157 0.9517721

#IPW
ps <- fitted(glmT)
wt <- data_exp$X/ps + (1 - data_exp$X)/(1 - ps)

summary(lm(Y ~ X+C1+C2, data = data_exp, weights = wt))$coef

# Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)  0.5736970 0.19133698  2.9983590 2.992716e-03
# X            0.9367955 0.13331067  7.0271609 2.061013e-11
# C1          -0.0252445 0.14051273 -0.1796598 8.575675e-01
# C2           0.3238391 0.06441795  5.0271564 9.592222e-07


fitCausal(dat = data_exp,
          formulas = forms_exp[-2],
          family = c(1,1,1)
)
# Error in convergence of optim
# log-likelihood:  -2695.41
# Z ~ C1 * C2
# est.   s.e. sandwich
# (intercept)  0.725 0.1183   0.1163
# C1          -0.070 0.1469   0.1426
# C2           0.663 0.0538   0.0525
# C1:C2        0.589 0.0664   0.0646
# residual s.e.:  0.937 0.041 0.0413
#
# Y ~ X * C1 * C2
# est.   s.e. sandwich
# (intercept)  0.9271 0.1897   0.1945
# X            0.4765 0.2339   0.2449
# C1          -0.3601 0.2795   0.2791
# C2           0.1905 0.0909   0.0994
# X:C1         0.3209 0.3221   0.3260
# X:C2         0.0521 0.1084   0.1197
# C1:C2        0.1778 0.1307   0.1332
# X:C1:C2     -0.0176 0.1473   0.1521
# residual s.e.:  0.987 0.0435 0.0461
#
# copula parameters:
#   cop ~ X + C1 + C2 + X:C1 + X:C2 + C1:C2 + X:C1:C2
# est.  s.e. sandwich
# (intercept)  0.4781 0.348    0.371
# C1           0.0149 0.574    0.627
# C2           0.1011 0.152    0.150
# X            0.5824 0.447    0.489
# C1:C2        0.1310 0.261    0.279
# C1:X        -0.6861 0.671    0.730
# C2:X        -0.1163 0.194    0.206
# C1:C2:X      0.1829 0.300    0.326
