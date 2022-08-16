// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(RcppClock)]]
#include <RcppArmadillo.h>
#include <RcppClock.h>

using namespace Rcpp;


//' @useDynLib ManyData
// static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
int g(int n) {
  if (n < 2) return(n);
  return(g(n-1) + g(n-2));
}

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// C++ function to compute density at points of Gaussian copula
// [[Rcpp::export]]
arma::vec dGcop(arma::mat const &x,
                arma::mat const &sigma,
                bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows;
  //d = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;

  //  double const constants = -(double)d/2.0 * log2pi;

  // arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  arma::vec const eig = arma::eig_sym(sigma);
  if (any(eig < 0)) {
    out = NA_REAL;
    // for (uword i = 0; i < n; i++) out(i) = nan;
    return out;
  }
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag()));

  // Rprintf("%3e\n", rootisum);

  for (uword i = 0; i < n; i++) {
    z = x.row(i);
    inplace_tri_mat_mult(z, rooti);
    //    Rprintf("%3e %3e\n", z(0), z(1));
    out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));
  }

  if (logd)
    return out;
  return exp(out);
}

// C++ function to compute density at points of Gaussian copula
//' @export
// [[Rcpp::export]]
arma::vec dGcop_sig(arma::mat const &x,
                    arma::cube const &sigma,
                    bool const logd = false) {

  Rcpp::Clock clock;

  using arma::uword;
  uword const n = x.n_rows;
  arma::vec out(n);
  // double const constants = -(double)d/2.0 * log2pi;

  for (uword i = 0; i < n; i++) {

    arma::vec eig = arma::eig_sym(sigma.slice(i));
    if (any(eig < 0)) {
      out(i) = NA_REAL;
      continue;
    }

    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma.slice(i))));
    double const rootisum = arma::sum(log(rooti.diag()));
    arma::rowvec z;

    // Rprintf("%3e\n", rootisum);

    z = x.row(i);

    inplace_tri_mat_mult(z, rooti);

    // Rprintf("%3e %3e\n", z(0), z(1));


    out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));

  }

  if (logd)
    return out;
  return exp(out);

}


// arma::vec dGcop_sig_clock(arma::mat const &x,
//                     arma::cube const &sigma,
//                     bool const logd = false) {
//
//   Rcpp::Clock clock;
//
//   using arma::uword;
//   uword const n = x.n_rows;
//   arma::vec out(n);
//   // double const constants = -(double)d/2.0 * log2pi;
//
//   for (uword i = 0; i < n; i++) {
//
//     clock.tick("eig_sym");
//     arma::vec eig = arma::eig_sym(sigma.slice(i));
//     if (any(eig < 0)) {
//       out(i) = NA_REAL;
//       continue;
//     }
//
//     clock.tock("eig_sym");
//
//     clock.tick("chol");
//     arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma.slice(i))));
//     clock.tock("chol");
//
//     clock.tick("sumlog");
//     double const rootisum = arma::sum(log(rooti.diag()));
//     clock.tock("sumlog");
//
//     arma::rowvec z;
//
//     // Rprintf("%3e\n", rootisum);
//
//     z = x.row(i);
//     clock.tick("inplace_tri_mat_mult");
//     inplace_tri_mat_mult(z, rooti);
//     clock.tock("inplace_tri_mat_mult");
//     // Rprintf("%3e %3e\n", z(0), z(1));
//
//     clock.tick("out_calc");
//     out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));
//     clock.tock("out_calc");
//   }
//
//   if (logd)
//     return out;
//   return exp(out);
//
//   clock.stop("clock");
// }
