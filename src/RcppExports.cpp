// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// qnormC2
arma::mat qnormC2(arma::mat& x);
RcppExport SEXP _ManyData_qnormC2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(qnormC2(x));
    return rcpp_result_gen;
END_RCPP
}
// qnormC3
arma::mat qnormC3(arma::mat& x);
RcppExport SEXP _ManyData_qnormC3(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(qnormC3(x));
    return rcpp_result_gen;
END_RCPP
}
// qnormC4
arma::mat qnormC4(arma::mat& x);
RcppExport SEXP _ManyData_qnormC4(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(qnormC4(x));
    return rcpp_result_gen;
END_RCPP
}
// createSig2
arma::cube createSig2(int ncv, int n, arma::mat par);
RcppExport SEXP _ManyData_createSig2(SEXP ncvSEXP, SEXP nSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ncv(ncvSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(createSig2(ncv, n, par));
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
    {"_ManyData_g", (DL_FUNC) &_ManyData_g, 1},
    {"_ManyData_dGcop", (DL_FUNC) &_ManyData_dGcop, 3},
    {"_ManyData_dGcop_sig", (DL_FUNC) &_ManyData_dGcop_sig, 3},
    {"_ManyData_expitC", (DL_FUNC) &_ManyData_expitC, 1},
    {"_ManyData_createSig", (DL_FUNC) &_ManyData_createSig, 3},
    {"_ManyData_qnormC", (DL_FUNC) &_ManyData_qnormC, 1},
    {"_ManyData_llC", (DL_FUNC) &_ManyData_llC, 5},
    {"_ManyData_qnormC2", (DL_FUNC) &_ManyData_qnormC2, 1},
    {"_ManyData_qnormC3", (DL_FUNC) &_ManyData_qnormC3, 1},
    {"_ManyData_qnormC4", (DL_FUNC) &_ManyData_qnormC4, 1},
    {"_ManyData_createSig2", (DL_FUNC) &_ManyData_createSig2, 3},
    {"_ManyData_univarDensC", (DL_FUNC) &_ManyData_univarDensC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ManyData(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}