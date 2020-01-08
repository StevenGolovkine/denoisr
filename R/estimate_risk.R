################################################################################
#            Functions for risk estimation using regularity                    #
################################################################################


#' Perform the estimation of the risk on a set of curves.
#'
#' @param curves List of (true) curves.
#' @param curves_estim List of estimated curves.
#' @param t0 The time point where to compute the risk
#' @return List with the mean and max residual squared error in t0.
#' @export
estimate_risk <- function(curves, curves_estim, t0 = 0.5) {
  risk <- estimateRisk(curves, curves_estim, t0)

  c("MeanRSE" = risk[1], "MaxRSE" = risk[2])
}

#' Perform the estimation of the risk on a set of curves along the points.
#'
#' @param curves List of (true) curves.
#' @param curves_estim List of estimated curves.
#'
#' @return List with the mean of and max the integrated residual squared error.
#' @export
estimate_risks <- function(curves, curves_estim) {
  risk <- estimateRiskCurves(curves, curves_estim)

  c("MeanIntRSE" = risk[1], "MaxIntRSE" = risk[2])
}
