################################################################################
#                     Functions for moments estimation                         #
################################################################################

#' Perform an estimation of the mean
#' 
#' This function performs the estimation of the mean of a set of curves.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U A vector of numerics, sampling points at which estimate the curves.
#' @param t0_list A vector of numerics, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0_list A vector of numerics, the number of neighbors of \eqn{t_0} to 
#'  consider. Should be set as \deqn{k0 = M* exp(-log(log(M))^2)}. We can set a 
#'  different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just 
#'  put a unique numeric.
#'  @return A list of with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points (equal to U)
#'   \item \strong{$x} The observed points.
#'  }
#' @export
mean_ll <- function(data, U = seq(0, 1, length.out = 101),
                    t0_list = 0.5, k0_list = 2,
                    grid = lseq(0.001, 0.1, length.out = 101),
                    nb_obs_minimal = 2, K = 'uniform'){
  
  data_smooth <- smooth_curves_mean(data, U = U, t0_list = t0_list,
                                    k0_list = k0_list, grid = grid,
                                    nb_obs_minimal = nb_obs_minimal, K = K)
  mu <- data_smooth$smooth %>% 
    purrr::map_dfc(~ .x$x) %>% 
    rowMeans(na.rm = TRUE)
  
  list(
    "parameter" = data_smooth$parameter,
    "mu" = mu
  )
}

#' Perform an estimation of the covariance
#' 
#' This function performs the estimation of the covariance of a set of curves.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U A vector of numerics, sampling points at which estimate the curves.
#' @param t0_list A vector of numerics, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0_list A vector of numerics, the number of neighbors of \eqn{t_0} to 
#'  consider. Should be set as \deqn{k0 = M* exp(-log(log(M))^2)}. We can set a 
#'  different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just 
#'  put a unique numeric.
#' @param centered Boolean, default=FALSE. Does the data have already been
#'  centered?
#'
#' @return A matrix of shape (length(U), length(U)) which is an estimation of
#'  the covariance matrix.
#' @export
covariance <- function(data, U = seq(0, 1, length.out = 101),
                       t0_list = 0.5, k0_list = 2,
                       t0_list_mean = seq(0.1, 0.9, by = 0.1), 
                       centered = FALSE, dense = FALSE){
  
  if(!centered){
    if(!dense){
      data_ <- list2cai(data)
      time_to_pred <- data_$time
      
      mean_curve <- mean_ll(data, U = time_to_pred, 
                            t0_list = t0_list_mean, k0_list = k0_list)
      data_$x <- data_$x - mean_curve
      data_unmean <- cai2list(data_$time, data_$x, data_$obs)
    } else {
      mean_curve <- mean_ll(data, U = data[[1]]$t, 
                            t0_list = t0_list_mean, k0_list = k0_list)
      data_unmean <- data %>% map(~ list(t = .x$t,
                                         x = .x$x - mean_curve))
    }
    data <- data_unmean
  }
  
  data_smooth <- smooth_curves(data, U = U, t0_list = t0_list,
                               k0_list = k0_list, reason = "covariance")
  data_smooth <- data_smooth$smooth
  
  cov_sum <- cov_count <- matrix(0, length(U), length(U))
  for(obs in 1:length(data_smooth)){
    obs_points <- which(!is.na(data_smooth[[obs]]$x))
    cov_count[obs_points, obs_points] <- cov_count[obs_points, obs_points] + 1
    cov_sum[obs_points, obs_points] <- cov_sum[obs_points, obs_points] + tcrossprod(data_smooth[[obs]]$x[obs_points])
  }
  cov_sum / cov_count
}

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
#' @param n The curve to remove from the estimation
#' @param t0_list A list of time at which the bandwidths have been estimated.
#'  (default = NULL)
#'  
#' @return A vector representing the LOO mean curve.
#' @export
#' @examples
#' df <- generate_fractional_brownian(N = 10, M = 300, H = 0.5, sigma = 0.05)
#' estimate_LOO_mean(df, b = 0.01, n = 1)
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
#' @examples
#' df <- generate_fractional_brownian(N = 10, M = 300, H = 0.5, sigma = 0.05)
#' estimate_mean(df, U = seq(0, 1, length.out = 101), b = 0.01)
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

#' Perform an estimation of the covariance
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
#' @examples
#' df <- generate_fractional_brownian(N = 10, M = 300, H = 0.5, sigma = 0.05)
#' estimate_covariance(df, U = seq(0, 1, length.out = 11), b = 0.01, h = 0.05)
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
