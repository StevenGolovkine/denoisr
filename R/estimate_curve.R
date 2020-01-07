################################################################################
#                 Functions that performs kernel smoothing                     #
################################################################################


#' Perform the smoothing of individuals curves.
#'
#' @importFrom magrittr %>%
#'
#' @param curve List of two vectors representing a curve:
#'               - $t Sampling points
#'               - $x Observed points
#' @param U Vector of points for the estimation
#' @param b Bandwith
#' @param t0_list List of time at which the bandwidths have been estimated.
#' @param kernel Kernel name to use
#'  - epanechnikov (default)
#'  - uniforma
#'  - beta
#'  - mBeta
#'  - locPoly (not used)
#' @useDynLib SmoothCurves
#'
#' @return List of two vectors representating the estimated curve:
#'               - $t Sampling points
#'               - $x Predicted points
#' @export
estimate_curve <- function(curve, U, b, t0_list = NULL,
                           kernel = "epanechnikov") {

  # Control of the bandwidth parameter
  if (length(b) == 1) {
    bandwidth <- rep(b, length(U))
  } else if (length(b) != length(U)) {
    idx <- U %>% purrr::map_dbl(~ order(abs(.x - t0_list))[1])
    bandwidth <- b[idx]
  } else if (length(b) == length(U)) {
    bandwidth <- b
  } else {
    stop("Issues with the bandwidth parameter.")
  }

  # Compute the estimation of the curve.
  if (kernel == "epanechnikov") {
    x_hat <- epaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else if (kernel == "uniform") {
    x_hat <- uniKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else if (kernel == "beta") {
    x_hat <- betaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else if (kernel == "mBeta") {
    x_hat <- modifiedBetaKernelSmoothingCurve(U, curve$t, curve$x, bandwidth)
  } else {
    print("Wrong kernel name")
    x_hat <- rep(0, length(U))
  }

  return(list(
    t = U,
    x = as.vector(x_hat)
  ))
}
