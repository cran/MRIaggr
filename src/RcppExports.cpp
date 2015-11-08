// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// filtrage2D_cpp
List filtrage2D_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, bool bilateral, bool na_rm);
RcppExport SEXP MRIaggr_filtrage2D_cpp(SEXP M_dataSEXP, SEXP M_operateurSEXP, SEXP index_dataSEXP, SEXP bilateralSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M_data(M_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M_operateur(M_operateurSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type index_data(index_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type bilateral(bilateralSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(filtrage2D_cpp(M_data, M_operateur, index_data, bilateral, na_rm));
    return __result;
END_RCPP
}
// filtrage3D_cpp
List filtrage3D_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, const NumericVector& Vec_operateur, const IntegerVector& p_operateur, const arma::mat& index_data, bool bilateral, bool na_rm);
RcppExport SEXP MRIaggr_filtrage3D_cpp(SEXP Vec_dataSEXP, SEXP p_dataSEXP, SEXP Vec_operateurSEXP, SEXP p_operateurSEXP, SEXP index_dataSEXP, SEXP bilateralSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type Vec_data(Vec_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_data(p_dataSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Vec_operateur(Vec_operateurSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_operateur(p_operateurSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type index_data(index_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type bilateral(bilateralSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(filtrage3D_cpp(Vec_data, p_data, Vec_operateur, p_operateur, index_data, bilateral, na_rm));
    return __result;
END_RCPP
}
// filtrage2Dmed_cpp
arma::mat filtrage2Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, bool na_rm);
RcppExport SEXP MRIaggr_filtrage2Dmed_cpp(SEXP M_dataSEXP, SEXP M_operateurSEXP, SEXP index_dataSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M_data(M_dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M_operateur(M_operateurSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type index_data(index_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(filtrage2Dmed_cpp(M_data, M_operateur, index_data, na_rm));
    return __result;
END_RCPP
}
// filtrage3Dmed_cpp
NumericVector filtrage3Dmed_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, const NumericVector& Vec_operateur, const IntegerVector& p_operateur, const arma::mat& index_data, bool na_rm);
RcppExport SEXP MRIaggr_filtrage3Dmed_cpp(SEXP Vec_dataSEXP, SEXP p_dataSEXP, SEXP Vec_operateurSEXP, SEXP p_operateurSEXP, SEXP index_dataSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type Vec_data(Vec_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_data(p_dataSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Vec_operateur(Vec_operateurSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_operateur(p_operateurSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type index_data(index_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(filtrage3Dmed_cpp(Vec_data, p_data, Vec_operateur, p_operateur, index_data, na_rm));
    return __result;
END_RCPP
}
// calcHemi_cpp
List calcHemi_cpp(const IntegerVector& coordsI, const IntegerVector& coordsJ, List ls_indexK, int n_num, const NumericVector& value, int n, double i_pos, double j_pos, double angle_pos, double penaltyNA, double sd_data, double p, bool symetrie);
RcppExport SEXP MRIaggr_calcHemi_cpp(SEXP coordsISEXP, SEXP coordsJSEXP, SEXP ls_indexKSEXP, SEXP n_numSEXP, SEXP valueSEXP, SEXP nSEXP, SEXP i_posSEXP, SEXP j_posSEXP, SEXP angle_posSEXP, SEXP penaltyNASEXP, SEXP sd_dataSEXP, SEXP pSEXP, SEXP symetrieSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector& >::type coordsI(coordsISEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type coordsJ(coordsJSEXP);
    Rcpp::traits::input_parameter< List >::type ls_indexK(ls_indexKSEXP);
    Rcpp::traits::input_parameter< int >::type n_num(n_numSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type i_pos(i_posSEXP);
    Rcpp::traits::input_parameter< double >::type j_pos(j_posSEXP);
    Rcpp::traits::input_parameter< double >::type angle_pos(angle_posSEXP);
    Rcpp::traits::input_parameter< double >::type penaltyNA(penaltyNASEXP);
    Rcpp::traits::input_parameter< double >::type sd_data(sd_dataSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type symetrie(symetrieSEXP);
    __result = Rcpp::wrap(calcHemi_cpp(coordsI, coordsJ, ls_indexK, n_num, value, n, i_pos, j_pos, angle_pos, penaltyNA, sd_data, p, symetrie));
    return __result;
END_RCPP
}
// calcContro_cpp
List calcContro_cpp(const arma::mat& contrast, const arma::mat& coords_px, const IntegerVector& index_k, const IntegerVector& index_k_contro, double d_lim, double lambda, int param_ref, double var_ref, bool type_moy, bool type_med, bool type_NN, bool verbose);
RcppExport SEXP MRIaggr_calcContro_cpp(SEXP contrastSEXP, SEXP coords_pxSEXP, SEXP index_kSEXP, SEXP index_k_controSEXP, SEXP d_limSEXP, SEXP lambdaSEXP, SEXP param_refSEXP, SEXP var_refSEXP, SEXP type_moySEXP, SEXP type_medSEXP, SEXP type_NNSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type contrast(contrastSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_px(coords_pxSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type index_k(index_kSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type index_k_contro(index_k_controSEXP);
    Rcpp::traits::input_parameter< double >::type d_lim(d_limSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type param_ref(param_refSEXP);
    Rcpp::traits::input_parameter< double >::type var_ref(var_refSEXP);
    Rcpp::traits::input_parameter< bool >::type type_moy(type_moySEXP);
    Rcpp::traits::input_parameter< bool >::type type_med(type_medSEXP);
    Rcpp::traits::input_parameter< bool >::type type_NN(type_NNSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(calcContro_cpp(contrast, coords_px, index_k, index_k_contro, d_lim, lambda, param_ref, var_ref, type_moy, type_med, type_NN, verbose));
    return __result;
END_RCPP
}
// calcMultiPotential_cpp
List calcMultiPotential_cpp(const S4& W_SR, const S4& W_LR, const NumericVector& sample, double& threshold, const arma::mat& coords, NumericVector& distance_ref, int& nbGroup_min, bool& multiV, double& neutre);
RcppExport SEXP MRIaggr_calcMultiPotential_cpp(SEXP W_SRSEXP, SEXP W_LRSEXP, SEXP sampleSEXP, SEXP thresholdSEXP, SEXP coordsSEXP, SEXP distance_refSEXP, SEXP nbGroup_minSEXP, SEXP multiVSEXP, SEXP neutreSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type W_SR(W_SRSEXP);
    Rcpp::traits::input_parameter< const S4& >::type W_LR(W_LRSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< double& >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type distance_ref(distance_refSEXP);
    Rcpp::traits::input_parameter< int& >::type nbGroup_min(nbGroup_minSEXP);
    Rcpp::traits::input_parameter< bool& >::type multiV(multiVSEXP);
    Rcpp::traits::input_parameter< double& >::type neutre(neutreSEXP);
    __result = Rcpp::wrap(calcMultiPotential_cpp(W_SR, W_LR, sample, threshold, coords, distance_ref, nbGroup_min, multiV, neutre));
    return __result;
END_RCPP
}
// calcPotts_cpp
List calcPotts_cpp(const S4& W_SR, const S4& W_LR, arma::mat sample, NumericVector rho, const arma::mat& coords, const IntegerVector& site_order, int iter_max, double cv_criterion, bool test_regional, NumericVector distance_ref, double threshold, double neutre, int nbGroup_min, bool multiV, bool last_vs_others, bool prior, bool type_reg, int verbose);
RcppExport SEXP MRIaggr_calcPotts_cpp(SEXP W_SRSEXP, SEXP W_LRSEXP, SEXP sampleSEXP, SEXP rhoSEXP, SEXP coordsSEXP, SEXP site_orderSEXP, SEXP iter_maxSEXP, SEXP cv_criterionSEXP, SEXP test_regionalSEXP, SEXP distance_refSEXP, SEXP thresholdSEXP, SEXP neutreSEXP, SEXP nbGroup_minSEXP, SEXP multiVSEXP, SEXP last_vs_othersSEXP, SEXP priorSEXP, SEXP type_regSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type W_SR(W_SRSEXP);
    Rcpp::traits::input_parameter< const S4& >::type W_LR(W_LRSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type site_order(site_orderSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< double >::type cv_criterion(cv_criterionSEXP);
    Rcpp::traits::input_parameter< bool >::type test_regional(test_regionalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_ref(distance_refSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type neutre(neutreSEXP);
    Rcpp::traits::input_parameter< int >::type nbGroup_min(nbGroup_minSEXP);
    Rcpp::traits::input_parameter< bool >::type multiV(multiVSEXP);
    Rcpp::traits::input_parameter< bool >::type last_vs_others(last_vs_othersSEXP);
    Rcpp::traits::input_parameter< bool >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< bool >::type type_reg(type_regSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(calcPotts_cpp(W_SR, W_LR, sample, rho, coords, site_order, iter_max, cv_criterion, test_regional, distance_ref, threshold, neutre, nbGroup_min, multiV, last_vs_others, prior, type_reg, verbose));
    return __result;
END_RCPP
}
// simulPotts_cpp
List simulPotts_cpp(const S4& W_SR, const S4& W_LR, arma::mat sample, const arma::mat& coords, const IntegerVector& site_order, NumericVector rho, NumericVector distance_ref, int iter_max, double cv_criterion, bool regional, bool verbose);
RcppExport SEXP MRIaggr_simulPotts_cpp(SEXP W_SRSEXP, SEXP W_LRSEXP, SEXP sampleSEXP, SEXP coordsSEXP, SEXP site_orderSEXP, SEXP rhoSEXP, SEXP distance_refSEXP, SEXP iter_maxSEXP, SEXP cv_criterionSEXP, SEXP regionalSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const S4& >::type W_SR(W_SRSEXP);
    Rcpp::traits::input_parameter< const S4& >::type W_LR(W_LRSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type site_order(site_orderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_ref(distance_refSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< double >::type cv_criterion(cv_criterionSEXP);
    Rcpp::traits::input_parameter< bool >::type regional(regionalSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(simulPotts_cpp(W_SR, W_LR, sample, coords, site_order, rho, distance_ref, iter_max, cv_criterion, regional, verbose));
    return __result;
END_RCPP
}
// simulPottsFast_cpp
arma::mat simulPottsFast_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const NumericVector& W_x, const IntegerVector& site_order, arma::mat sample, double rho, int n, const int p, int iter_nb);
RcppExport SEXP MRIaggr_simulPottsFast_cpp(SEXP W_iSEXP, SEXP W_pSEXP, SEXP W_xSEXP, SEXP site_orderSEXP, SEXP sampleSEXP, SEXP rhoSEXP, SEXP nSEXP, SEXP pSEXP, SEXP iter_nbSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_i(W_iSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_p(W_pSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type W_x(W_xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type site_order(site_orderSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type iter_nb(iter_nbSEXP);
    __result = Rcpp::wrap(simulPottsFast_cpp(W_i, W_p, W_x, site_order, sample, rho, n, p, iter_nb));
    return __result;
END_RCPP
}
// calcGroupsCoords_cpp
List calcGroupsCoords_cpp(const arma::mat& coords_NNA, const IntegerVector& index_NNA, const arma::mat& Neighborhood, IntegerVector& coords_max, int& max_groups, bool& verbose);
RcppExport SEXP MRIaggr_calcGroupsCoords_cpp(SEXP coords_NNASEXP, SEXP index_NNASEXP, SEXP NeighborhoodSEXP, SEXP coords_maxSEXP, SEXP max_groupsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_NNA(coords_NNASEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type index_NNA(index_NNASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Neighborhood(NeighborhoodSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type coords_max(coords_maxSEXP);
    Rcpp::traits::input_parameter< int& >::type max_groups(max_groupsSEXP);
    Rcpp::traits::input_parameter< bool& >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(calcGroupsCoords_cpp(coords_NNA, index_NNA, Neighborhood, coords_max, max_groups, verbose));
    return __result;
END_RCPP
}
// calcGroupsW_cpp
List calcGroupsW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& subset, int& max_groups);
RcppExport SEXP MRIaggr_calcGroupsW_cpp(SEXP W_iSEXP, SEXP W_pSEXP, SEXP subsetSEXP, SEXP max_groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_i(W_iSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_p(W_pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< int& >::type max_groups(max_groupsSEXP);
    __result = Rcpp::wrap(calcGroupsW_cpp(W_i, W_p, subset, max_groups));
    return __result;
END_RCPP
}
// calcRadius_cpp
List calcRadius_cpp(const arma::mat& coords, const NumericVector& sample, double& threshold, const LogicalVector& subset_bary, bool& verbose);
RcppExport SEXP MRIaggr_calcRadius_cpp(SEXP coordsSEXP, SEXP sampleSEXP, SEXP thresholdSEXP, SEXP subset_barySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< double& >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type subset_bary(subset_barySEXP);
    Rcpp::traits::input_parameter< bool& >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(calcRadius_cpp(coords, sample, threshold, subset_bary, verbose));
    return __result;
END_RCPP
}
// calcBlockW_cpp
List calcBlockW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& site_order, const NumericVector& dist_center, double dist_max, bool& verbose);
RcppExport SEXP MRIaggr_calcBlockW_cpp(SEXP W_iSEXP, SEXP W_pSEXP, SEXP site_orderSEXP, SEXP dist_centerSEXP, SEXP dist_maxSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_i(W_iSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type W_p(W_pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type site_order(site_orderSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type dist_center(dist_centerSEXP);
    Rcpp::traits::input_parameter< double >::type dist_max(dist_maxSEXP);
    Rcpp::traits::input_parameter< bool& >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(calcBlockW_cpp(W_i, W_p, site_order, dist_center, dist_max, verbose));
    return __result;
END_RCPP
}
