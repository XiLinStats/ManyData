//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends(RcppClock)]]
#include <RcppArmadillo.h>
#include <RcppClock.h>
#include <RcppDist.h>
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
// //' @export
// // [[Rcpp::export]]
// arma::vec llC_clock(DataFrame dat,
//              arma::mat mm,
//              arma::mat beta,
//              arma::vec phi,
//              arma::vec inCop
//                 ){
//
//   Rcpp::Clock clock;
//
//
//   arma::mat log_den(dat.nrow(),dat.ncol());
//   arma::mat dat_u(dat.nrow(),dat.ncol());
//
//   clock.tick("eta");
//   arma::mat eta = mm * beta;
//   clock.tock("eta");
//   // Rcout << "The value of eta : " << eta << "\n";
//
//   int n = dat.ncol();
//
//   clock.tick("univarDens");
//
//   for (int i = 0 ; i <n ; ++i){
//     arma::vec eta_i = eta.col(i);
//
//     // NumericVector check = dat[i];
//
//     // Rcout << "The value of dat[i] : " << check << "\n";
//     List L = univarDensC(dat[i],as<NumericVector>(wrap(eta_i)),phi[i]);
//     log_den.col(i) = as<arma::vec>(L["ld"]);
//     dat_u.col(i) = as<arma::vec>(L["u"]);
//   }
//
//   clock.tock("univarDens");
//
//   clock.tick("create par matrix");
//
//   int nv = phi.size();
//   int ncv = inCop.size();
//   int n_comb = Rf_choose(ncv,2);
//   // Rcout << "The value of n_comb : " << n_comb << "\n";
//
//   arma::mat par(dat.nrow(),n_comb);
//
//   for (int i = 0 ; i < n_comb ; ++i){
//     NumericVector out =  pmin(pmax(2 *expitC(as<NumericVector>(wrap(eta.col(i+nv)))) - 1, -1 + 1e-10), 1 - 1e-10);
//     par.col(i) = as<arma::vec>(out);
//   }
//
//   clock.tock("create par matrix");
//
//   // Create sigma cube
//
//   clock.tick("create sigma cube");
//   arma::cube sigma = createSig(ncv,dat.nrow(),par);
//   clock.tock("create sigma cube");
//
//   clock.tick("qnorm");
//   arma::mat qx = qnormC(dat_u);
//   clock.tock("qnorm");
//
//   clock.tick("dGcop_sig");
//   arma::vec cop = dGcop_sig(qx, sigma, true);
//   clock.tock("dGcop_sig");
//
//   clock.tick("summing");
//   arma::vec out = cop + arma::sum(log_den,1);
//   clock.tock("summing");
//
//   clock.stop("clock");
//
//   return out;
// }
//
// //
//' @export
// [[Rcpp::export]]
arma::vec llC(DataFrame dat,
              arma::mat mm,
              arma::mat beta,
              arma::vec phi,
              arma::vec inCop
){




  arma::mat log_den(dat.nrow(),dat.ncol());
  arma::mat dat_u(dat.nrow(),dat.ncol());


  arma::mat eta = mm * beta;

  // Rcout << "The value of eta : " << eta << "\n";

  int n = dat.ncol();



  for (int i = 0 ; i <n ; ++i){
    arma::vec eta_i = eta.col(i);

    // NumericVector check = dat[i];

    // Rcout << "The value of dat[i] : " << check << "\n";
    List L = univarDensC(dat[i],as<NumericVector>(wrap(eta_i)),phi[i]);
    log_den.col(i) = as<arma::vec>(L["ld"]);
    dat_u.col(i) = as<arma::vec>(L["u"]);
  }

  int nv = phi.size();
  int ncv = inCop.size();
  int n_comb = Rf_choose(ncv,2);
  // Rcout << "The value of n_comb : " << n_comb << "\n";

  arma::mat par(dat.nrow(),n_comb);

  for (int i = 0 ; i < n_comb ; ++i){
    NumericVector out =  pmin(pmax(2 *expitC(as<NumericVector>(wrap(eta.col(i+nv)))) - 1, -1 + 1e-10), 1 - 1e-10);
    par.col(i) = as<arma::vec>(out);
  }

  arma::cube sigma = createSig(ncv,dat.nrow(),par);
  arma::mat qx = qnormC(dat_u);
  arma::vec cop = dGcop_sig(qx, sigma, true);
  arma::vec out = cop + arma::sum(log_den,1);



  return out;
}



//' @export
// [[Rcpp::export]]
double postll_C(
                arma::rowvec current_val,

                const arma::vec inCop,
                const DataFrame dat_exp,
                const arma::vec theta_exp,
                const arma::mat mm_exp,
                const List mask_exp,

                const DataFrame dat_obs,
                const arma::vec theta_obs,
                const arma::mat mm_obs,
                const List mask_obs,

                const arma::vec p_mu,
                const arma::mat p_sigma,

                const double eta
){

  // initial values
  arma::vec theta_exp2 = theta_exp;
  arma::vec theta_obs2 = theta_obs;
  theta_exp2(2) = theta_obs2(2) = current_val(0);
  theta_exp2(3) = theta_obs2(3) = current_val(1);
  theta_exp2(4) = theta_obs2(4) = current_val(2);
  theta_exp2(5) = theta_obs2(5) = current_val(3);


  // fill beta matrix and phi vector with theta values

  arma::mat beta_exp = mask_exp["beta_m"];

  int c = 0;

  for(auto& val : beta_exp){
    if (val == 1){
      val = theta_exp2[c];
      c += 1;
    }

  }

  arma::vec phi_exp = mask_exp["phi_m"];

  for(auto& val : phi_exp){
    if (val == 1){
      val = theta_exp2[c];
      c += 1;
    }

  }

  arma::mat beta_obs = mask_obs["beta_m"];

  int k = 0;

  for(auto& val : beta_obs){
    if (val == 1){
      val = theta_obs2[k];
      k += 1;
    }

  }

  arma::vec phi_obs = mask_obs["phi_m"];

  for(auto& val : phi_obs){
    if (val == 1){
      val = theta_obs2[k];
      k += 1;
    }
//
  }

  double post = sum(llC(dat_exp,mm_exp,beta_exp, phi_exp,inCop)) + eta * sum(llC(dat_obs, mm_obs,beta_obs, phi_obs,inCop)) + arma::sum(dmvnorm(current_val, p_mu, p_sigma, true));

  // const arma::vec p_mu = {0,0};
  // const arma::mat p_sigma = {{2,0},{0,2}};
  // Rcout << "The value of p_sigma : " << p_sigma << "\n";


  return post;



}




// [[Rcpp::export]]
arma::mat MCMCloop_C(
    const arma::uword n_iter,

    const arma::rowvec init_val,
    const arma::mat sigma,

    const arma::vec inCop,
    const DataFrame dat_exp,
    arma::vec theta_exp,
    const arma::mat mm_exp,
    const List mask_exp,

    const DataFrame dat_obs,
    arma::vec theta_obs,
    const arma::mat mm_obs,
    const List mask_obs,

    const arma::vec p_mu,
    const arma::mat p_sigma,

    const double eta
){


  arma::mat chain(n_iter,4);
  arma::vec u(n_iter,arma::fill::randu);

  chain.row(0) = init_val;

  double prev_postll = postll_C(init_val,inCop,
                                 dat_exp,theta_exp,mm_exp,mask_exp,
                                 dat_obs,theta_obs,mm_obs,mask_obs,
                                 p_mu,p_sigma,eta);

  arma::rowvec prev_var = init_val;

  for (arma::uword i = 1 ; i < n_iter ; ++i){
    arma::rowvec new_var = rmvnorm(1, prev_var.t(),sigma);

    double curr_postll = postll_C(new_var,inCop,
                                  dat_exp,theta_exp,mm_exp,mask_exp,
                                  dat_obs,theta_obs,mm_obs,mask_obs,
                                  p_mu,p_sigma,eta);

    double alpha = exp( curr_postll - prev_postll);

    if (u(i) < alpha){
      chain.row(i) = new_var;
      prev_var = new_var;
      prev_postll = curr_postll;
    }else{chain.row(i) = prev_var;
    }



  }

  return(chain);


}



