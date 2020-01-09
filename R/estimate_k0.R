################################################################################
#         Functions for k0 parameter estimation using regularity               #
################################################################################

#' Perform the estimation of the pilot \eqn{k_0}
#'
#' This function performs the estimation of the pilot \eqn{k_0} as used in 
#' \cite{add ref}.
#' 
#' @family estimate \eqn{k_0}
#' 
#' @param M Numeric, mean number of sampling points per curve
#'
#' @return Numeric, the pilot \eqn{k_0}
#' @export
#' @examples 
#' estimate_k0_pilot(200)
estimate_k0_pilot <- function(M) {
  trunc((M + 1) * exp(-sqrt(log(M + 1)))) + 1
}

#' Perform the estimation of the oracle \eqn{k_0}
#' 
#' This function performs the estimation of the oracle \eqn{k_0} as used in
#' \cite{add ref}.
#' 
#' @family estimate \eqn{k_0}
#'
#' @param M Numeric, mean number of sampling points per curve
#' @param H Numeric, estimation of \eqn{H_0}
#'
#' @return Numeric, the oracle \eqn{k_0}
#' @export
#' @examples 
#' estimate_k0_oracle(200, 0.5)
estimate_k0_oracle <- function(M, H) {
  trunc((M + 1)**(H / (2 + H))) + 1
}
