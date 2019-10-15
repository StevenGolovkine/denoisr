######################################################################################
#                  Functions for risk estimation using regularity                    #
######################################################################################
library(tidyverse)

Rcpp::sourceCpp('./src/estimate_risk.cpp')


#' Perform the estimation of the risk on a set of curves. 
#' 
#' @param curves List of (true) curves.
#' @param curves_estim List of estimated curves.
#' @param t0 The time point where to compute the risk
#' @return List with the mean and max residual squared error in t0.
estimate.risk <- function(curves, curves_estim, t0 = 0.5){
  
  # Estimate the risk
  risk <- estimateRisk(curves, curves_estim, t0)
  
  #return(c('MeanRSE' = risk[1], 'MaxRSE' = risk[2]))
  return(risk)
}

#' Perform the estimation of the risk on a set of curves along the points.
#' 
#' @param curves List of (true) curves.
#' @param curves_estim List of estimated curves.
#' 
#' @return List with the mean of mean and max of max residual squared error.
estimate.risks <- function(curves, curves_estim){
  
  # Estimate the risk
  risk <- estimateRiskCurves(curves, curves_estim)
  
  return(c('MeanRSE' = risk[1], 'MaxRSE' = risk[2]))
}