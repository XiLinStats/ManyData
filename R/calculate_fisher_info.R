#' Calculate the sum of log-likelihood
#'
#' @param masks A masks object from masks().
#' @param theta A theta vector.
#' @param mm A model matrix from causl::mm().
#' @param dat Input data to calculate log-likelihood on.
#' @return The sum of log-likelihood.
#' @examples
#'
#'

ll_sum <- function(masks, theta, mm, dat){

  theta2 <- copy(theta)
  msks2 <- copy(masks)
  np <- sum(msks2$beta_m > 0)
  msks2$beta_m[masks$beta_m > 0] <- theta2[seq_len(np)]
  msks2$phi_m[masks$phi_m > 0] <- theta2[-seq_len(np)]
  ll_i <- causl:::ll(dat = dat[,c(1,5)], mm = mm, beta = msks2$beta_m, phi = msks2$phi_m, inCop = c(1,2),
                     fam_cop = 1, family = list(1,1), link = NULL, par2 = NULL,
                     useC = TRUE)
  return(sum(ll_i))
}


#' Numerically approximate the fisher information matrix using finite difference method
#'
#' @param masks A masks object from masks().
#' @param theta A theta vector.
#' @param mm A model matrix from causl::mm().
#' @param dat Input data to calculate log-likelihood on.
#' @param delta The infinitesimal delta used in finite difference method.
#' @return The sum of log-likelihood.
#' @examples
#'
#'
ApproxFI <- function(masks, theta, mm, dat, delta = 1e-4){


  d_11<- (ll_sum(masks, Add_Delta(theta, 6, 2*delta), mm ,dat) - 2 * ll_sum(masks, Add_Delta(theta, 6, delta), mm ,dat) + ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_12 <- (ll_sum(masks, Add_Delta(theta, c(6,8), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(6), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(8), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_22 <- (ll_sum(masks, Add_Delta(theta, 8, 2*delta), mm ,dat) -  2*ll_sum(masks, Add_Delta(theta, c(8), delta) , mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  FI<- -matrix(data = c(d_11, d_12, d_12, d_22),nrow =2)

  return(FI)
}


ApproxFI_all <- function(masks, theta, mm, dat, delta = 1e-4){

  d_11<- (ll_sum(masks, Add_Delta(theta, 6, 2*delta), mm ,dat) - 2 * ll_sum(masks, Add_Delta(theta, 6, delta), mm ,dat) + ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_12 <- (ll_sum(masks, Add_Delta(theta, c(6,8), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(6), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(8), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_13 <- (ll_sum(masks, Add_Delta(theta, c(6,5), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(6), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(5), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_14 <- (ll_sum(masks, Add_Delta(theta, c(6,7), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(6), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(7), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_22 <- (ll_sum(masks, Add_Delta(theta, 8, 2*delta), mm ,dat) -  2*ll_sum(masks, Add_Delta(theta, c(8), delta) , mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_23 <- (ll_sum(masks, Add_Delta(theta, c(8,5), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(8), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(5), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_24 <- (ll_sum(masks, Add_Delta(theta, c(8,7), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(8), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(7), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_33<- (ll_sum(masks, Add_Delta(theta, 5, 2*delta), mm ,dat) - 2 * ll_sum(masks, Add_Delta(theta, 5, delta), mm ,dat) + ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_34 <- (ll_sum(masks, Add_Delta(theta, c(5,7), delta), mm ,dat) -  ll_sum(masks, Add_Delta(theta, c(5), delta), mm ,dat) - ll_sum(masks, Add_Delta(theta, c(7), delta), mm ,dat) +  ll_sum(masks, theta, mm ,dat))/(delta^2)

  d_44<- (ll_sum(masks, Add_Delta(theta, 7, 2*delta), mm ,dat) - 2 * ll_sum(masks, Add_Delta(theta, 7, delta), mm ,dat) + ll_sum(masks, theta, mm ,dat))/(delta^2)

  FI<- -matrix(data = c(d_11, d_12, d_13, d_14,
                        d_12, d_22, d_23, d_24,
                        d_13, d_23, d_33, d_34,
                        d_14, d_24, d_34, d_44),nrow =4)

  return(FI)
}

Add_Delta <- function(theta, list, delta){
  theta[list] <- theta[list] + delta
  return(theta)
}
