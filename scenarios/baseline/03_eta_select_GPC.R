# install.packages("parallel")
library(parallel)


# Parallelisation ---------------------------------------------------------

parallel::detectCores()
# [1] 12
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()



# Start the loop ----------------------------------------------------------

# True parameters
theta_true <- c(0.2,0.8)

# Number of bootstrap samples
n_boot <- 500

# List of eta
eta_list <- seq(0.6,0.6,0.1)

# Function to evaluate

set.seed(111)

func <- function(eta) {

  cover <- matrix(0,1,2)
  data_exp_b <- data_exp[sample(nrow(data_exp), nrow(data_exp),replace = T), ]
  data_obs_b <- data_obs[sample(nrow(data_obs), nrow(data_obs),replace = T), ]

  theta_sim_raw <- run_MH_MCMC(c(1,0.5), 2000, dat_obs=data_obs_b, dat_exp = data_exp_b, eta = eta)
  theta_sim_rtnd <- theta_sim_raw[-seq(1:50),][seq(1, nrow(theta_sim_raw) - 50, by =1), ]
  #
  l1 <- quantile(theta_sim_rtnd[,1],probs = 0.025)
  u1 <-  quantile(theta_sim_rtnd[,1],probs = 0.975)

  l2 <- quantile(theta_sim_rtnd[,2],probs = 0.025)
  u2 <-  quantile(theta_sim_rtnd[,2],probs = 0.975)

  cover[1,] <- c((l1 < theta_true[1] & u1 > theta_true[1]) *1, (l2<theta_true[2] & u2 > theta_true[2]))

  return(cover)
}

# Initialise
cover_prob <- matrix(0,length(eta_list),3)
start_time <- Sys.time()

# Loop

for (i in 1 : length(eta_list)) {

  eta <- eta_list[[i]]
  print(paste0("Current Eta:",eta))

  final <- foreach(
    i = 1:n_boot,
    .combine = 'rbind',
    .packages = c('mvtnorm','data.table')
  ) %dopar% {
    x<-func(eta)
  }
  cover_prob[i,] <- c(eta,colMeans(final))
}

Sys.time() -  start_time










































