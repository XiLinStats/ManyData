
# Read in simulated dataset -----------------------------------------------

data_exp <- readRDS(file.path(data_dir,"/data_exp.rds"))
data_obs <- readRDS(file.path(data_dir,"/data_obs.rds"))

# Formulae
forms_obs <- list(Z ~ C1,
                  list(X ~ Z*C1,C1~ 1),
                  Y ~ X*C1,
                  ~ X*C1 #copula
)


forms_exp <- list(Z ~ C1,
                  list(X ~ 1,C1~ 1),
                  Y ~ X*C1,
                  ~ X*C1 #copula
)


# MLE from fitCausl ------------------------------------------------------

# Exp data only
fit_exp<-fitCausal(dat = data_exp,
                   formulas = forms_exp[-2],
                   family = c(1,1,1)
)
MLE_exp<- fit_exp$pars$Y$beta

MLE_exp[c(2,3)]<- fit_exp$pars$Y$beta[c(3,2)]
names(MLE_exp)[c(2,3)] <- c("C1","X")

# Obs data only
fit_obs <- fitCausal(dat = data_obs,
                   formulas = forms_obs[-2],
                   family = c(1,1,1)
)

MLE_obs<- fit_obs$pars$Y$beta

MLE_obs[c(2,3)]<- fit_obs$pars$Y$beta[c(3,2)]
names(MLE_obs)[c(2,3)] <- c("C1","X")

# Combine both data
# MLE from obs + exp combined
data_comb <- rbind(data_obs,data_exp)
fit_comb<-fitCausal(dat = data_comb,
                    formulas = forms_exp[-2],
                    family = c(1,1,1)
)

MLE_comb<- fit_comb$pars$Y$beta
MLE_comb[c(2,3)]<- fit_comb$pars$Y$beta[c(3,2)]
names(MLE_comb)[c(2,3)] <- c("C1","X")
