//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends(RcppClock)]]
#include <RcppArmadillo.h>
#include <RcppClock.h>
#include "shared.h"
using namespace Rcpp;
//
// // [[Rcpp::export]]
// arma::mat qnormC2( arma::mat& x) {
//
//   arma::mat out(x.n_rows,x.n_cols);
//   out = x;
//
//   for(auto& val : out){
//
//
//       val = R::qnorm(val,0,1,true,false);
//
//   }
//
//   return(out);
//
// }
//
//
//
// // [[Rcpp::export]]
// arma::mat qnormC3( arma::mat& x) {
//
//   arma::mat out(x.n_rows,x.n_cols);
//
//
//   for(auto& val : x){
//
//
//     val = R::qnorm(val,0,1,true,false);
//
//   }
//
//   //
//   //   int n = x.n_cols;
//   //   int n_row = x.n_rows;
//   //
//   //   for (int j = 0; j < n ; ++j){
//   //
//   //     for (int i = 0; i < n_row ; ++i){
//   //
//   //       out(i,j) = R::qnorm(x(i,j),0,1,true,false);
//   //     }
//   //
//   //   }
//
//
//   return(x);
//
// }
//
// // [[Rcpp::export]]
// arma::mat qnormC4( arma::mat& x) {
//
//   arma::mat out(x.n_rows,x.n_cols);
//
//
//     int n = x.n_cols;
//     int n_row = x.n_rows;
//
//     for (int j = 0; j < x.n_cols ; ++j){
//
//       for (int i = 0; i < n_row ; ++i){
//
//         out(i,j) = R::qnorm(x(i,j),0,1,true,false);
//       }
//
//     }
//
//
//   return(x);
//
// }
//
//
//
// // [[Rcpp::export]]
// arma::cube createSig2(int ncv, int n , arma::mat par) {
//
//   arma::cube sigma(ncv,ncv,n,arma::fill::ones);
//
//   for (int i = 0 ; i <n ; ++i){
//
//     int c = 0;
//
//     for (int j = 1;  j <ncv ; ++j){
//
//
//         for (int k = 0;  k <j ; ++k){
//
//           sigma(k,j,i) = sigma(j,k,i) =par(i,c);
//
//           c = c+1;
//
//
//
//       }
//     }
//
//
//   }
//   return sigma;
// }

// [[Rcpp::export]]
arma::mat expit2(arma::mat x){
  arma::mat out(x.n_rows,x.n_cols);
  out = x;
  for(auto& val : out){
  val = fmin(fmax(2 *(exp(val)/(1+exp(val))) - 1, -1 + 1e-10), 1 - 1e-10);
  }
  return out;
}


NumericVector expitC(NumericVector x) {
  NumericVector res = exp(x)/(1+exp(x));
  return res;
}

// [[Rcpp::export]]
arma::mat expit3(arma::mat x){

  arma::mat res(x.n_rows,x.n_cols);

  for (int i = 0 ; i < x.n_cols ; ++i){
    NumericVector out =  pmin(pmax(2 *expitC(as<NumericVector>(wrap(x.col(i)))) - 1, -1 + 1e-10), 1 - 1e-10);
    res.col(i) = as<arma::vec>(out);
  }
  return(res);

}


