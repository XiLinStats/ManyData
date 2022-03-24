//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends(RcppClock)]]
#include <RcppArmadillo.h>
#include <RcppClock.h>
#include "shared.h"
using namespace Rcpp;

//' @useDynLib ManyData
//' @import RcppArmadillo causl



// [[Rcpp::export]]
NumericVector expitC(NumericVector x) {
  NumericVector res = exp(x)/(1+exp(x));
  return res;
}


// [[Rcpp::export]]
arma::cube createSig(int ncv, int n , arma::mat par) {

  arma::cube sigma(ncv,ncv,n);

  for (int i = 0 ; i <n ; ++i){

    int c = 0;

    for (int j = 0;  j <ncv ; ++j){

      sigma(j,j,i) = 1;

      if (j > 0){

        for (int k = 0;  k <j ; ++k){

          sigma(k,j,i) = sigma(j,k,i) =par(i,c);

          c = c+1;

        }

      }
    }


  }
  return sigma;
}


// [[Rcpp::export]]
arma::mat qnormC(const arma::mat& x) {

  arma::mat out(x.n_rows,x.n_cols);
  out = x;

  for(auto& val : out){
    val = R::qnorm(val,0,1,true,false);
   }

  return(out);

}

//
//' @export
// [[Rcpp::export]]
arma::vec llC(DataFrame dat,
             arma::mat mm,
             arma::mat beta,
             arma::vec phi,
             arma::vec inCop
                ){

  Rcpp::Clock clock;


  arma::mat log_den(dat.nrow(),dat.ncol());
  arma::mat dat_u(dat.nrow(),dat.ncol());

  clock.tick("eta");
  arma::mat eta = mm * beta;
  clock.tock("eta");
  // Rcout << "The value of eta : " << eta << "\n";

  int n = dat.ncol();

  clock.tick("univarDens");

  for (int i = 0 ; i <n ; ++i){
    arma::vec eta_i = eta.col(i);

    // NumericVector check = dat[i];

    // Rcout << "The value of dat[i] : " << check << "\n";
    List L = univarDensC(dat[i],as<NumericVector>(wrap(eta_i)),phi[i]);
    log_den.col(i) = as<arma::vec>(L["ld"]);
    dat_u.col(i) = as<arma::vec>(L["u"]);
  }

  clock.tock("univarDens");

  clock.tick("create par matrix");

  int nv = phi.size();
  int ncv = inCop.size();
  int n_comb = Rf_choose(ncv,2);
  // Rcout << "The value of n_comb : " << n_comb << "\n";

  arma::mat par(dat.nrow(),n_comb);

  for (int i = 0 ; i < n_comb ; ++i){
    NumericVector out =  pmin(pmax(2 *expitC(as<NumericVector>(wrap(eta.col(i+nv)))) - 1, -1 + 1e-10), 1 - 1e-10);
    par.col(i) = as<arma::vec>(out);
  }

  clock.tock("create par matrix");

  // Create sigma cube

  clock.tick("create sigma cube");
  arma::cube sigma = createSig(ncv,dat.nrow(),par);
  clock.tock("create sigma cube");

  clock.tick("qnorm");
  arma::mat qx = qnormC(dat_u);
  clock.tock("qnorm");

  clock.tick("dGcop_sig");
  arma::vec cop = dGcop_sig(qx, sigma, true);
  clock.tock("dGcop_sig");

  clock.tick("summing");
  arma::vec out = cop + arma::sum(log_den,1);
  clock.tock("summing");

  clock.stop("clock");

  return out;
}
