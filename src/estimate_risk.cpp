// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_risk.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector estimateRisk(
    const List & curves, // Curves list ($x and $t)
    const List & curves_estim, // Estimated curves list ($x and $t)
    const double & t0 // In which point estimate the risk.
){
  // Define risk list
  List risk(2);
  
  // Get the number of curves
  arma::uword N = curves.length();
  
  List mycurve = curves[0];
  List mycurve_estim = curves_estim[0];
  
  arma::uvec idx;
  
  arma::vec squared_diff(N);
  for(arma::uword n=0; n<N; n++){
    
    mycurve = curves[n];
    arma::vec x = mycurve["x"];
    arma::vec t = mycurve["t"];
    
    mycurve_estim = curves_estim[n];
    arma::vec x_estim = mycurve_estim["x"];
    arma::vec t_estim = mycurve_estim["t"];
    
    double idx = arma::max(find(t <= t0));
    
    double real = x(idx);
    double pred = x_estim(idx);
    
    squared_diff(n) = std::pow(pred - real, 2);
  }
  
  risk(0) = arma::mean(squared_diff);
  risk(1) = arma::max(squared_diff);
  
  Rcpp::NumericVector tmp = Rcpp::wrap(squared_diff);
  tmp.attr("dim") = R_NilValue;
  
  return (tmp);
}

// [[Rcpp::export]]
List estimateRiskCurve(
    const List & curve, // Curve ($x and $t)
    const List & curve_estim // Estimated curve ($x and $t)
){
  List risk(2);
  
  arma::vec x = curve["x"];
  arma::vec x_estim = curve_estim["x"];
  
  arma::vec squared_diff = pow(x - x_estim, 2);
  risk(0) = arma::mean(squared_diff);
  risk(1) = arma::max(squared_diff);
  
  return (risk);
}

// [[Rcpp::export]]
List estimateRiskCurves(
    const List & curves, // Curves list ($x and $t)
    const List & curves_estim // Estimated curves list ($x and $t)
){
  List risk(2);
  
  // Get the number of curves
  arma::uword N = curves.length();
  
  List mycurve = curves[0];
  List mycurve_estim = curves_estim[0];
  
  arma::vec mean_risk(N);
  arma::vec max_risk(N);
  List risk_curve(2);
  for(arma::uword n; n < N; n++){
    
    mycurve = curves[n];
    mycurve_estim = curves_estim[n];

    risk_curve = estimateRiskCurve(mycurve, mycurve_estim);
    
    mean_risk(n) = risk_curve(0);
    max_risk(n) = risk_curve(1);
  }
  
  //risk(0) = arma::mean(mean_risk);
  //risk(1) = arma::max(max_risk);
  risk(0) = mean_risk;
  risk(1) = max_risk;
  
  Rcpp::NumericVector tmp = Rcpp::wrap(mean_risk);
  tmp.attr("dim") = R_NilValue;
  
  Rcpp::NumericVector tmp2 = Rcpp::wrap(max_risk);
  tmp2.attr("dim") = R_NilValue;
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("Mean") = tmp,
    Rcpp::Named("Max") = tmp2);
  
  return (result);
}