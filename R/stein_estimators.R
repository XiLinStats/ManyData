##' Shrinkage estimators
##'
##' Functions to implement shrinkage estimators of Green and
##' Strawderman (1991) and Rosenman et al. (2020)
##'
##' @param theta_e,theta_o estimates from the experimental and observational populations
##' @param sigma2_e estimates of the standard error for experimental coefficients
##'
##' @export
gs <- function(theta_e, theta_o, sigma2_e=rep(1,N)) {

  if (!is.matrix(theta_e)) theta_e <- matrix(theta_e, nrow=1)
  if (!is.matrix(theta_o)) theta_o <- matrix(theta_o, nrow=1)

  N <- nrow(theta_e)
  d <- ncol(theta_e)
  if (nrow(theta_o) != N) stop("Different number of values")
  if (ncol(theta_o) != d) stop("Different dimensions")
  if (length(sigma2_e) != N) stop("Variance estimates do not match number of observations")


  if (d <= 2) return(theta_e)

  scl <- pmax(0, (1 - sigma2_e*(d-2)/(rowSums((theta_e-theta_o)^2))))

  theta_o + scl*(theta_e - theta_o)
}

##' @describeIn gs Multivariate extension delta1
##' @param Sigma_e estimated covariance for \code{theta_e}
##' @param a coefficient for shrinkage
##' @export
gsar1 <- function(theta_e, theta_o, Sigma_e=diag(p), a=p-2) {

  # if (!is.matrix(theta_e)) theta_e <- matrix(theta_e, nrow=1)
  # if (!is.matrix(theta_o)) theta_o <- matrix(theta_o, nrow=1)

  # N <- nrow(theta_e)
  p <- length(theta_e)
  # if (nrow(theta_o) != N) stop("Different number of values")
  if (length(theta_o) != p) stop("Different dimensions")
  if (nrow(Sigma_e) != ncol(Sigma_e)) stop("Covariance must be square")
  if (nrow(Sigma_e) != p) stop("Variance estimates do not match number of variables")

  if (p <= 2) return(theta_e)

  scl <- pmax(0, (1 - a/c( t(theta_e - theta_o) %*% solve(Sigma_e) %*% (theta_e - theta_o) )))

  theta_o + scl*(theta_e - theta_o)
}

##' @describeIn gs Multivariate extension delta2
##' @export
gsar2 <- function(theta_e, theta_o, Sigma_e=diag(p), a=p-2) {

  # if (!is.matrix(theta_e)) theta_e <- matrix(theta_e, nrow=1)
  # if (!is.matrix(theta_o)) theta_o <- matrix(theta_o, nrow=1)

  # N <- nrow(theta_e)
  p <- length(theta_e)
  # if (nrow(theta_o) != N) stop("Different number of values")
  if (length(theta_o) != p) stop("Different dimensions")
  if (nrow(Sigma_e) != ncol(Sigma_e)) stop("Covariance must be square")
  if (nrow(Sigma_e) != p) stop("Variance estimates do not match number of variables")

  if (p <= 2) return(theta_e)

  iSigma <- solve(Sigma_e)
  iSigma2 <- iSigma %*% iSigma

  scl <- pmax(0, (1 - a*iSigma/c( t(theta_e - theta_o) %*% iSigma2 %*% (theta_e - theta_o) )))
  dim(scl) <- c(p,p)

  theta_o + scl %*% (theta_e - theta_o)
}



##' @describeIn gs Rosenman et al. (2020) kappa1 estimate
##' @export
rbob_kappa1 <- function (theta_e, theta_o, sigma2_e=rep(1,p)) {

  ## get length of vectors
  p <- length(theta_e)
  if (length(theta_o) != p) stop("Different number of strata in experimental and observational data")
  if (length(sigma2_e) != p) stop("Number of variance estimates does not match number of strata")

  # kappa_1+
  lambda_1 <- sum(sigma2_e)/sum((theta_e - theta_o)^2)
  kappa1_plus <- theta_o +
    max(1 - lambda_1, 0)*(theta_e - theta_o)

  return(kappa1_plus)
}

##' @describeIn gs Rosenman et al. (2020) kappa2 estimate
##' @export
rbob_kappa2 <- function (theta_e, theta_o, sigma2_e=rep(1,p)) {

  ## get length of vectors
  p <- length(theta_e)
  if (length(theta_o) != p) stop("Different number of strata in experimental and observational data")
  if (length(sigma2_e) != p) stop("Number of variance estimates does not match number of strata")

  # # kappa_1+
  # lambda_1 <- sum(sigma2_e)/sum((theta_e - theta_o)^2)
  # kappa1_plus <- theta_o +
  #   max(1 - lambda_1, 0)*(theta_e - theta_o)

  # kappa_2+
  lambda_2 <- pmin(pmax(sum(sigma2_e)^2/
                          sum(sigma2_e^2*(theta_e - theta_o)^2)*diag(sigma2_e), 0), 1)
  kappa2_plus <- c(theta_o +
                     pmax(diag(p) - sum(sigma2_e)^2/
                            sum(sigma2_e^2*(theta_e - theta_o)^2)*diag(sigma2_e), 0) %*%
                     (theta_e - theta_o))
  kappa2_plus <- c(theta_o +
                     pmax(diag(p) - lambda_2, 0) %*%
                     (theta_e - theta_o))


  return(kappa2_plus)
}
