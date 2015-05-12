// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// disclapglm_linkfun
NumericVector disclapglm_linkfun(NumericVector mu);
RcppExport SEXP disclapmix_disclapglm_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        NumericVector __result = disclapglm_linkfun(mu);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// disclapglm_linkinv
NumericVector disclapglm_linkinv(NumericVector eta);
RcppExport SEXP disclapmix_disclapglm_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP );
        NumericVector __result = disclapglm_linkinv(eta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// disclapglm_mu_eta
NumericVector disclapglm_mu_eta(NumericVector eta);
RcppExport SEXP disclapmix_disclapglm_mu_eta(SEXP etaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP );
        NumericVector __result = disclapglm_mu_eta(eta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// disclapglm_varfunc
NumericVector disclapglm_varfunc(NumericVector mu);
RcppExport SEXP disclapmix_disclapglm_varfunc(SEXP muSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        NumericVector __result = disclapglm_varfunc(mu);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// disclapglm_loglikeh
double disclapglm_loglikeh(double mu, double y);
RcppExport SEXP disclapmix_disclapglm_loglikeh(SEXP muSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type mu(muSEXP );
        Rcpp::traits::input_parameter< double >::type y(ySEXP );
        double __result = disclapglm_loglikeh(mu, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// disclapglm_deviance
double disclapglm_deviance(NumericVector y, NumericVector mu, NumericVector wt);
RcppExport SEXP disclapmix_disclapglm_deviance(SEXP ySEXP, SEXP muSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP );
        double __result = disclapglm_deviance(y, mu, wt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities_sum
NumericVector rcpp_calculate_haplotype_probabilities_sum(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP disclapmix_rcpp_calculate_haplotype_probabilities_sum(SEXP allele_rangeSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type allele_range(allele_rangeSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        NumericVector __result = rcpp_calculate_haplotype_probabilities_sum(allele_range, y, p, tau);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_match_quantities
NumericVector rcpp_match_quantities(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP disclapmix_rcpp_match_quantities(SEXP allele_rangeSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type allele_range(allele_rangeSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        NumericVector __result = rcpp_match_quantities(allele_range, y, p, tau);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_create_design_matrix
IntegerMatrix rcpp_create_design_matrix(IntegerMatrix x, int clusters);
RcppExport SEXP disclapmix_rcpp_create_design_matrix(SEXP xSEXP, SEXP clustersSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type clusters(clustersSEXP );
        IntegerMatrix __result = rcpp_create_design_matrix(x, clusters);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_create_new_weight_vector
NumericVector rcpp_create_new_weight_vector(NumericMatrix vic, int loci);
RcppExport SEXP disclapmix_rcpp_create_new_weight_vector(SEXP vicSEXP, SEXP lociSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type vic(vicSEXP );
        Rcpp::traits::input_parameter< int >::type loci(lociSEXP );
        NumericVector __result = rcpp_create_new_weight_vector(vic, loci);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_create_response_vector
IntegerVector rcpp_create_response_vector(IntegerMatrix x, IntegerMatrix y);
RcppExport SEXP disclapmix_rcpp_create_response_vector(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        IntegerVector __result = rcpp_create_response_vector(x, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_calculate_wic
NumericMatrix rcpp_calculate_wic(IntegerMatrix x, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP disclapmix_rcpp_calculate_wic(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        NumericMatrix __result = rcpp_calculate_wic(x, y, p, tau);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_calculate_vic
NumericMatrix rcpp_calculate_vic(NumericMatrix wic);
RcppExport SEXP disclapmix_rcpp_calculate_vic(SEXP wicSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type wic(wicSEXP );
        NumericMatrix __result = rcpp_calculate_vic(wic);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities
NumericVector rcpp_calculate_haplotype_probabilities(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau);
RcppExport SEXP disclapmix_rcpp_calculate_haplotype_probabilities(SEXP new_dataSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type new_data(new_dataSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        NumericVector __result = rcpp_calculate_haplotype_probabilities(new_data, y, p, tau);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_calculate_haplotype_probabilities_increase_at_alleles
NumericVector rcpp_calculate_haplotype_probabilities_increase_at_alleles(IntegerMatrix new_data, IntegerMatrix y, NumericMatrix p, NumericVector tau, NumericVector constants, IntegerVector increase_at_alleles);
RcppExport SEXP disclapmix_rcpp_calculate_haplotype_probabilities_increase_at_alleles(SEXP new_dataSEXP, SEXP ySEXP, SEXP pSEXP, SEXP tauSEXP, SEXP constantsSEXP, SEXP increase_at_allelesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type new_data(new_dataSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type constants(constantsSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type increase_at_alleles(increase_at_allelesSEXP );
        NumericVector __result = rcpp_calculate_haplotype_probabilities_increase_at_alleles(new_data, y, p, tau, constants, increase_at_alleles);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_simulate
IntegerMatrix rcpp_simulate(int nsim, IntegerMatrix y, NumericVector tau_cumsum, NumericMatrix disclap_parameters);
RcppExport SEXP disclapmix_rcpp_simulate(SEXP nsimSEXP, SEXP ySEXP, SEXP tau_cumsumSEXP, SEXP disclap_parametersSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau_cumsum(tau_cumsumSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type disclap_parameters(disclap_parametersSEXP );
        IntegerMatrix __result = rcpp_simulate(nsim, y, tau_cumsum, disclap_parameters);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_find_haplotype_in_matrix
int rcpp_find_haplotype_in_matrix(const IntegerMatrix subpop, const IntegerVector h);
RcppExport SEXP disclapmix_rcpp_find_haplotype_in_matrix(SEXP subpopSEXP, SEXP hSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const IntegerMatrix >::type subpop(subpopSEXP );
        Rcpp::traits::input_parameter< const IntegerVector >::type h(hSEXP );
        int __result = rcpp_find_haplotype_in_matrix(subpop, h);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
