######################################################################################
#                     Functions that performs kernel smoothing                       #
######################################################################################
library(tidyverse)

Rcpp::sourceCpp('./src/estimate_curve.cpp')

#' Perform the smoothing of individuals curves
#'
#' @param curve List of two vectors representing a curve:
#'               - $t Sampling points
#'               - $x Observed points
#' @param U Vector of points for the estimation
#' @param b Bandwith 
#' @param kernel Kernel name to use, default='epanechnikov'
#' @useDynLib CovarianceEstimate
#' @return List of two vectors representaing the estimated curve:
#'               - $t Sampling points
#'               - $x Predicted points
estimate.curve <- function(curve, U, b, kernel='epanechnikov', degree = 0){
  
  # Compute the estimation of the curve.
  if(kernel == 'epanechnikov'){
    x_hat <- kernelSmoothingCurve(U, curve$t, curve$x, b)
  } else if(kernel == 'beta'){
    x_hat <- betaKernelSmoothingCurve(U, curve$t, curve$x, b)
  } else if(kernel == 'mBeta'){
    x_hat <- modifiedBetaKernelSmoothingCurve(U, curve$t, curve$x, b)
  #} else if(kernel == 'locPoly'){
    #require(KernSmooth)
    #smooth <- locpoly(curve$t, curve$x, bandwidth = b, 
    #                 gridsize = 2*length(U), degree = degree)
    #x_hat <- smooth$y[U == smooth$x]
  } else{
    print("Wrong kernel name")
    x_hat <- rep(0, length(U))
  }
  
  return(list(t = U,
              x = as.vector(x_hat)))
}
