################################################################################
#         Functions for k0 parameter estimation using regularity               #
################################################################################


#' Perform the estimation of the pilot k0.
#'
#' @param M Mean number of samping points per curve.
#'
#' @return The pilot k0.
estimate_k0_pilot <- function(M) {
  k0_pilot <- trunc((M + 1) * exp(-sqrt(log(M + 1)))) + 1
  return(k0_pilot)
}

#' Perform the estimation of the oracle k0.
#'
#' @param M Mean number of sampling points per curve.
#' @param H Estimation of the regularity of the curves.
#'
#' @return The oracle k0.
estimate_k0_oracle <- function(M, H) {
  k0_oracle <- trunc((M + 1)**(H / (2 + H))) + 1
  return(k0_oracle)
}
