// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// expitC
NumericVector expitC(NumericVector x);
RcppExport SEXP _ManyData_expitC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expitC(x));
    return rcpp_result_gen;
END_RCPP
}
// createSig
arma::cube createSig(int ncv, int n, arma::mat par);
RcppExport SEXP _ManyData_createSig(SEXP ncvSEXP, SEXP nSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ncv(ncvSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(createSig(ncv, n, par));
    return rcpp_result_gen;
END_RCPP
}
// qnormC
arma::mat qnormC(const arma::mat& x);
RcppExport SEXP _ManyData_qnormC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(qnormC(x));
    return rcpp_result_gen;
END_RCPP
}
// llC
arma::vec llC(DataFrame dat, arma::mat mm, arma::mat beta, arma::vec phi, arma::vec inCop);
RcppExport SEXP _ManyData_llC(SEXP datSEXP, SEXP mmSEXP, SEXP betaSEXP, SEXP phiSEXP, SEXP inCopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type inCop(inCopSEXP);
    rcpp_result_gen = Rcpp::wrap(llC(dat, mm, beta, phi, inCop));
    return rcpp_result_gen;
END_RCPP
}
// postll_C_all
double postll_C_all(arma::rowvec current_val, const arma::vec inCop, const DataFrame dat_exp, const arma::vec theta_exp, const arma::mat mm_exp, const List mask_exp, const DataFrame dat_obs, const arma::vec theta_obs, const arma::mat mm_obs, const List mask_obs, const arma::vec p_mu, const arma::mat p_sigma, const double eta);
RcppExport SEXP _ManyData_postll_C_all(SEXP current_valSEXP, SEXP inCopSEXP, SEXP dat_expSEXP, SEXP theta_expSEXP, SEXP mm_expSEXP, SEXP mask_expSEXP, SEXP dat_obsSEXP, SEXP theta_obsSEXP, SEXP mm_obsSEXP, SEXP mask_obsSEXP, SEXP p_muSEXP, SEXP p_sigmaSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type current_val(current_valSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type inCop(inCopSEXP);
    Rcpp::traits::input_parameter< const DataFrame >::type dat_exp(dat_expSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta_exp(theta_expSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type mm_exp(mm_expSEXP);
    Rcpp::traits::input_parameter< const List >::type mask_exp(mask_expSEXP);
    Rcpp::traits::input_parameter< const DataFrame >::type dat_obs(dat_obsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta_obs(theta_obsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type mm_obs(mm_obsSEXP);
    Rcpp::traits::input_parameter< const List >::type mask_obs(mask_obsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type p_mu(p_muSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type p_sigma(p_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(postll_C_all(current_val, inCop, dat_exp, theta_exp, mm_exp, mask_exp, dat_obs, theta_obs, mm_obs, mask_obs, p_mu, p_sigma, eta));
    return rcpp_result_gen;
END_RCPP
}
// MCMCloop_C
arma::mat MCMCloop_C(const arma::uword n_iter, const arma::rowvec init_val, const arma::vec sigma, const arma::vec inCop, const DataFrame dat_exp, arma::vec theta_exp, const arma::mat mm_exp, const List mask_exp, const DataFrame dat_obs, arma::vec theta_obs, const arma::mat mm_obs, const List mask_obs, const arma::vec p_mu, const arma::mat p_sigma, const double eta);
RcppExport SEXP _ManyData_MCMCloop_C(SEXP n_iterSEXP, SEXP init_valSEXP, SEXP sigmaSEXP, SEXP inCopSEXP, SEXP dat_expSEXP, SEXP theta_expSEXP, SEXP mm_expSEXP, SEXP mask_expSEXP, SEXP dat_obsSEXP, SEXP theta_obsSEXP, SEXP mm_obsSEXP, SEXP mask_obsSEXP, SEXP p_muSEXP, SEXP p_sigmaSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type inCop(inCopSEXP);
    Rcpp::traits::input_parameter< const DataFrame >::type dat_exp(dat_expSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_exp(theta_expSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type mm_exp(mm_expSEXP);
    Rcpp::traits::input_parameter< const List >::type mask_exp(mask_expSEXP);
    Rcpp::traits::input_parameter< const DataFrame >::type dat_obs(dat_obsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_obs(theta_obsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type mm_obs(mm_obsSEXP);
    Rcpp::traits::input_parameter< const List >::type mask_obs(mask_obsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type p_mu(p_muSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type p_sigma(p_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMCloop_C(n_iter, init_val, sigma, inCop, dat_exp, theta_exp, mm_exp, mask_exp, dat_obs, theta_obs, mm_obs, mask_obs, p_mu, p_sigma, eta));
    return rcpp_result_gen;
END_RCPP
}
// g
int g(int n);
RcppExport SEXP _ManyData_g(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(g(n));
    return rcpp_result_gen;
END_RCPP
}
// dGcop
arma::vec dGcop(arma::mat const& x, arma::mat const& sigma, bool const logd);
RcppExport SEXP _ManyData_dGcop(SEXP xSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dGcop(x, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// dGcop_sig
arma::vec dGcop_sig(arma::mat const& x, arma::cube const& sigma, bool const logd);
RcppExport SEXP _ManyData_dGcop_sig(SEXP xSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::cube const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dGcop_sig(x, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// univarDensC
List univarDensC(NumericVector x, NumericVector eta, double phi, String link, int family, int df);
RcppExport SEXP _ManyData_univarDensC(SEXP xSEXP, SEXP etaSEXP, SEXP phiSEXP, SEXP linkSEXP, SEXP familySEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< int >::type family(familySEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(univarDensC(x, eta, phi, link, family, df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ManyData_expitC", (DL_FUNC) &_ManyData_expitC, 1},
    {"_ManyData_createSig", (DL_FUNC) &_ManyData_createSig, 3},
    {"_ManyData_qnormC", (DL_FUNC) &_ManyData_qnormC, 1},
    {"_ManyData_llC", (DL_FUNC) &_ManyData_llC, 5},
    {"_ManyData_postll_C_all", (DL_FUNC) &_ManyData_postll_C_all, 13},
    {"_ManyData_MCMCloop_C", (DL_FUNC) &_ManyData_MCMCloop_C, 15},
    {"_ManyData_g", (DL_FUNC) &_ManyData_g, 1},
    {"_ManyData_dGcop", (DL_FUNC) &_ManyData_dGcop, 3},
    {"_ManyData_dGcop_sig", (DL_FUNC) &_ManyData_dGcop_sig, 3},
    {"_ManyData_univarDensC", (DL_FUNC) &_ManyData_univarDensC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ManyData(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
