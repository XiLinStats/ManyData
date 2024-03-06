##' Comparator methods
##'
##' Functions to implement comparator methods.
##' Currently, Oberst et al. (2022)
##'
##' @param theta_e,theta_o estimates from the experimental and observational populations
##' @param var_e,var_o estimates of the standard error for the estimates from the experimental and observational populations
##'
##' @export
obserst <- function(theta_e, theta_o, var_e, var_o){

  lambda <- (var_e)/((theta_e - theta_o)^2 + var_e + var_o)
  theta_lambda <- lambda * theta_o + (1 - lambda) * theta_e

  return(list(theta_oberst = theta_lambda, lambda = lambda))
}
