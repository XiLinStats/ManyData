using namespace Rcpp;
Rcpp::NumericVector timesTwo(NumericVector x);

List univarDensC(NumericVector x, NumericVector eta, double phi , String link = "missing", int family = 1 ,
                 int df = 50);

arma::vec dGcop(arma::mat const &x,
                arma::mat const &sigma,
                bool const logd = false);

arma::vec dGcop_sig(arma::mat const &x,
                    arma::cube const &sigma,
                    bool const logd = false);
