#' SmoothCurves: A package for the smoothing of a set of curves.
#' 
#' The package provides functions in order to smooth a set of curves. The idea is
#' to perform non-parametric kernel regression of each of curves using the 
#' Nadaraya-Watson estimator. In order to use this estimator, we need to select
#' a bandwidth. A common way to select this bandwidth is to used cross-validation.
#' However, we are going to use the large number of curves to estimate it more 
#' efficiently. The idea is to estimate the underlying regularity of the curves
#' (assuming they have the same). Then, the estimation of the bandwidth follows.
#' 
#' @docType package
#' @name SmoothCurves
NULL