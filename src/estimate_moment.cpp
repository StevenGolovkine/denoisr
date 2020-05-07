// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_mean.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "estimate_curve.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec LOOmean(
  const List & curves, // Curves list ($x and $t)
  const arma::vec & U, // Estimation points
  const List & b, // Smoothing bandwidths
  const arma::uword & n // Curve to remove from the estimation
){
  
  // Get parameters
  arma::uword N = curves.length(); // Number of curves
  arma::uword L = U.n_elem; // Number of sampling points
  
  List mycurve = curves[0];
  
  // Define output
  arma::mat res(L, N - 1);
  res.fill(0);
  
  arma::uword cpt = 0;
  for(arma::uword m=0; m<N; m++){
    if(m != n){
      mycurve = curves[m];
      arma::vec X = mycurve["x"];
      arma::vec T = mycurve["t"];
      
      arma::vec bandwidth = b[m];

      res.col(cpt) = epaKernelSmoothingCurve(U, T, X, bandwidth);
      
      cpt++;
    }
  }
  
  return mean(res, 1);
}


// [[Rcpp::export]]
arma::mat covariance(
    const List & curves, // Curves list ($x and $t)
    const arma::vec & sampling_points, // Estimation points
    const List & b, // Smoothing bandwidth 
    const List & h // Smoothing bandwidth
){
  
  // Get parameters
  arma::uword N = curves.length(); // Number of curves
  arma::vec U = sampling_points; // Number of sampling points
  List mycurve = curves[0];
  
  // Define output
  arma::mat res(U.n_elem, U.n_elem);
  
  arma::mat smooths_U(U.n_elem, N);
  arma::mat smooths_V(U.n_elem, N);
  for(arma::uword n=0; n<N; n++){
    mycurve = curves[n];
    arma::vec T = mycurve["t"];
    arma::vec X = mycurve["x"];
    
    arma::vec mean_curve = LOOmean(curves, T, b, n);
    arma::vec smooth_curve = epaKernelSmoothingCurve(T, T, X, b[n]);
    arma::vec smooth_curve_unmean = smooth_curve - mean_curve;
    
    smooths_U.col(n) = epaKernelSmoothingCurve(U, T, smooth_curve_unmean, h[n]);
    smooths_V.col(n) = epaKernelSmoothingCurve(U, T, smooth_curve_unmean, h[n]);
  }

  res = cov(smooths_U.t(), smooths_V.t(), 1); // Normalisation using N
  return res;
}
