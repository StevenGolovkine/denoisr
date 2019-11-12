// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_risk.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List estimateRisk(
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
  
  return (risk);
}

// [[Rcpp::export]]
double estimateRiskCurve(
    const List & curve, // Curve ($x and $t)
    const List & curve_estim // Estimated curve ($x and $t)
){
  arma::mat risk;
  
  arma::vec x = curve["x"];
  arma::vec t = curve["t"];
  arma::vec x_estim = curve_estim["x"];
  
  arma::vec squared_diff = pow(x - x_estim, 2);
  risk = arma::trapz(t, squared_diff); // Numerical integration
  
  return (risk(0, 0));
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
  
  arma::vec risk_int(N);
  for(arma::uword n=0; n < N; n++){
    
    mycurve = curves[n];
    mycurve_estim = curves_estim[n];

    risk_int(n) = estimateRiskCurve(mycurve, mycurve_estim);
  }
  
  risk(0) = arma::mean(risk_int);
  risk(1) = arma::max(risk_int);
  
  return (risk);
}