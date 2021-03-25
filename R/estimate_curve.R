################################################################################
#                 Functions that performs kernel smoothing                     #
################################################################################


#' Perform the smoothing of individual curve.
#'
#' This function performs the smoothing of a curve using the Nadaraya-Watson 
#' estimator given a particular kernel.
#' 
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp
#' 
#' @param curve A list, with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U A vector of numerics, sampling points at which estimate the curve. 
#' @param b Numeric or vector of numerics, estimation of the bandwidth. If one 
#'  is provided, we use a unique bandwidth for the curve. However, if a vector 
#'  is given, the bandwidth changes depending on the sampling points. 
#' @param t0_list A vector of numerics, times at which the bandwidths have been 
#'  estimated. Only used if the parameter \code{b} is a vector.
#' @param kernel Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item uniform
#'   \item beta
#'   \item mBeta 
#'  }
#' @param n_obs_min Integer, minimum number of observation for the smoothing
#' @useDynLib denoisr
#'
#' @return A list, with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  } 
#' @export
#' @examples
#' df_piece <- generate_piecewise_fractional_brownian(N = 1, M = 300, 
#'                                                    H = c(0.2, 0.5, 0.8), 
#'                                                    sigma = 0.05)
#' estimate_curve(df_piece[[1]],
#'                U = seq(0, 1, length.out = 200),
#'                b = c(0.2, 0.5, 0.8),
#'                t0_list = c(0.16, 0.5, 0.83))
estimate_curve <- function(curve, U, b, t0_list = NULL,
                           kernel = "epanechnikov", n_obs_min = 1) {
  if (length(b) == 1) {
    bandwidth <- rep(b, length(U))
  } else if ((length(b) != length(U)) & !is.null(t0_list)) {
    #idx <- U %>% purrr::map_dbl(~ order(abs(.x - t0_list))[1])
    #bandwidth <- b[idx]
    bandwidth <- approx(t0_list, b, xout = U,
                        yleft = b[1], yright = b[length(b)])$y
  } else if (length(b) == length(U)) {
    bandwidth <- b
  } else {
    stop("Issues with the bandwidth parameter.")
  }

  if (kernel == "epanechnikov") {
    x_hat <- epaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "uniform") {
    x_hat <- uniKernelSmoothingCurve(U, curve$t, curve$x, bandwidth, n_obs_min)
  } else if (kernel == "beta") {
    x_hat <- betaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else if (kernel == "mBeta") {
    x_hat <- modifiedBetaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else {
    print("Wrong kernel name")
    x_hat <- rep(0, length(U))
  }

  list(t = U, x = as.vector(x_hat))
}
