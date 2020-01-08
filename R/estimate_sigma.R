################################################################################
#       Functions for sigma parameter estimation using regularity              #
################################################################################


#' Perform the estimation of sigma
#'
#' @param data List of curves to estimate by kernel regression.
#' @return An estimation of sigma.
#' @export
estimate_sigma <- function(data) {
  estimateSigma(data)
}
