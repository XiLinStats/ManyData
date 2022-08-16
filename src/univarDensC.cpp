#include <Rcpp.h>
using namespace Rcpp;

//' @useDynLib ManyData
//' @import Rcpp
//' @export
// [[Rcpp::export]]
List univarDensC(NumericVector x, NumericVector eta, double phi , String link = "missing", int family = 1 ,
                 int df = 50){

  int n = x.size();

  NumericVector ld(n);
  NumericVector u(n);
  NumericVector mu(n);
  double sigma = sqrt(phi);

  // Gaussian;

  if (family == 1){

    if (link == "missing" || link == "identity"){
      mu = eta;
    }else if(link == "inverse"){
      mu = pow(eta,-1);
    }else if (link == "log") {
      mu = exp(eta);
      // Rcout << "The value of mu : " << mu << "\n";
    }else {
      stop("Not a valid link function for Gaussian distribution");
    }

    for (int i = 0; i < n; ++i){

      ld[i] = R::dnorm(x[i],mu[i],sigma,true);
      u[i] = R::pnorm(x[i],mu[i],sigma , true, false);
      // Rcout << "The value of out[i] : " << out[i] << "\n";
    }


    // t-distribution;
  } else if (family == 2){

    if (link == "missing" || link == "identity"){
      mu = eta;
    }else if(link == "inverse"){
      mu = pow(eta,-1);
    }else if (link == "log") {
      mu = exp(eta);
      // Rcout << "The value of mu : " << mu << "\n";
    }else {
      stop("Not a valid link function for t-distribution");
    }

    for (int i = 0; i < n; ++i){

      ld[i] = R::dt((x[i] - mu[i])/sigma, df = df, true) -
        log(sigma);
      u[i] = R::pt((x[i] - mu[i])/sigma, df = df, true, false);
      // Rcout << "The value of out[i] : " << out[i] << "\n";
    }


    //  Gamma distribution;
  } else if (family == 3){
    // log is the default link
    if (link == "missing" || link == "log"){
      mu = exp(eta);
    }else if(link == "inverse"){
      mu = pow(eta,-1);
    }else if (link == "identity") {
      mu = eta;

    }else {
      stop("Not a valid link function for gamma distribution");
    }

    for (int i = 0; i < n; ++i){

      ld[i] = R::dgamma(x[i], 1/phi,  phi * mu[i], true);
      u[i] = R::pgamma(x[i],  1/phi,  phi * mu[i], true, false);
      // Rcout << "The value of out[i] : " << out[i] << "\n";
    }

  }else {
    stop("Only Gaussian, t, gamma and Bernoulli distributions are allowed");
  }
  List L = List::create(Named("ld") = ld,Named("u")= u);

  return L;
}



