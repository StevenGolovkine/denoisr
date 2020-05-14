################################################################################
#                     Functions for moments estimation                         #
################################################################################

#' Perform a leave-one-out mean curve
#' 
#' This function performs the computation of a leave-one-out mean curve.
#' 
#' @param curves A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param b A list of bandwidth
#' @param n The curves to remove from the estimation
#' @param t0_list A list of time at which the bandwidths have been estimated.
#'  (default = NULL)
#'  
#' @return A vector representing the LOO mean curve.
#' @export
estimate_LOO_mean <- function(curves, b, n, t0_list = NULL){
  
  U <- curves[[n]][['t']]
  
  if (length(b) == 1) {
    bandwidth <- purrr::rerun(length(curves), rep(b, length(U)))
  } else if ((length(b) == length(curves)) & (length(b[[1]]) == 1)) {
    bandwidth <- purrr::map(b, ~ rep(.x, length(U)))
  } else if ((length(b) == length(curves)) & length(b[[1]] > 1)) {
    idx <- U %>% purrr::map_dbl(~ order(abs(.x - t0_list))[1])
    bandwidth <- b %>% purrr::map(~ .x[idx])
  } else {
    stop("Issues with the bandwidth parameter.")
  }
  
  as.vector(LOOmean(curves, U, bandwidth, n - 1))
}

#' Perform an estimation of the leave-one-out mean curve
#' 
#' This function performs the estimation of the leave-one-out mean curve of a
#' set of curves.
#' 
#' @param curves A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param U A vector of numerics, sampling points at which estimate the curve.
#' @param b A list of bandwidth
#' @param t0_list A list of time at which the bandwidths have been estimated.
#'  (default = NULL)
#'  
#' @return A vector representing the mean curve.
#' @export
estimate_mean <- function(curves, U, b, t0_list = NULL){
  
  if ((length(b[[1]]) > 1) & is.null(t0_list)){
    stop(paste("If there are multiple bandwidths for each cuves,",
               "t0_list can not be NULL."))
  }
  
  if (length(b) == 1) {
    bandwidth <- purrr::rerun(length(curves), rep(b, length(U)))
  } else if ((length(b) == length(curves)) & (length(b[[1]]) == 1)) {
    bandwidth <- purrr::map(b, ~ rep(.x, length(U)))
  } else if ((length(b) == length(curves)) & length(b[[1]] > 1)) {
    idx <- U %>% purrr::map_dbl(~ order(abs(.x - t0_list))[1])
    bandwidth <- b %>% purrr::map(~ .x[idx])
  } else {
    stop("Issues with the bandwidth parameter.")
  }
  
  as.vector(mean_cpp(curves, U, bandwidth))
} 

#' Perform an estimation of the covariancce
#' 
#' This function performs the estimation of the covariance of a set of curves.
#'
#' @param curves A list, where each element represents a real curve. Each curve
#'  have to be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param U A vector of numerics, sampling points at which estimate the curve.
#' @param b A list of bandwidth for the first smoothing
#' @param h A list of bandwidth for the second smoothing
#' @param t0_list A list of time at which the bandwidths have been estimated.
#'  (default = NULL)
#'  
#' @return A matrix of shape (length(U), length(U)) which is an estimation of
#'  the covariance matrix.
#' @export
estimate_covariance <- function(curves, U, b, h, t0_list = NULL){
  
  
  if (length(b) == 1) {
    band_b <- curves %>% purrr::map(~ rep(b, length(.x$t)))
  } else if ((length(b) == length(curves)) & (length(b[[1]]) == 1)) {
    band_b <- purrr::map2(b, curves, ~ rep(.x, length(.y$t)))
  } else if ((length(b) == length(curves)) & length(b[[1]] > 1)) {
    
    idx <- curves %>% purrr::map(~ purrr::map_dbl(.x$t, ~ order(abs(.x - t0_list))[1]))
    band_b <- purrr::map2(b, idx, ~ .x[.y])
  } else {
    stop("Issues with the bandwidth b parameter.")
  }
  
  if (length(h) == 1) {
    band_h <- purrr::rerun(length(curves), rep(h, length(U)))
  } else if ((length(h) == length(curves)) & (length(h[[1]]) == 1)) {
    band_h <- purrr::map(h, ~ rep(.x, length(U)))
  } else if ((length(h) == length(curves)) & length(h[[1]] > 1)) {
    idx <- U %>% purrr::map_dbl(~ order(abs(.x - t0_list))[1])
    band_h <- h %>% purrr::map(~ .x[idx])
  } else {
    stop("Issues with the bandwidth h parameter.")
  }
  
  LOOmean_list <- purrr::map(1:length(curves), 
                             ~ estimate_LOO_mean(curves, b, .x, t0_list))
  covariance_cpp(curves, LOOmean_list, U, band_b, band_h)
}
