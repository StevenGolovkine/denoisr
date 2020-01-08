################################################################################
#         Functions for H0 parameter estimation using regularity               #
################################################################################


#' Perform the estimation of H0.
#'
#' @importFrom magrittr %>%
#'
#' @param data List of curves to estimate by kernel regression.
#' @param t0 The starting time for the estimation of H0. We consider the 8k0 - 7
#'  nearest points of t0 for the estimation of H0 when sigma is unknown.
#' @param k0 For the computation of the gap between the different observations.
#'  Should be set as k0 = M / max(8, log(M)).
#' @param sigma True value of sigma
#'  If null, change estimate.
#'
#' @return An estimation of H0.
estimate_H0 <- function(data, t0 = 0, k0 = 2, sigma = NULL) {
  theta <- function(v, k, idx) (v[idx + 2 * k - 1] - v[idx + k])**2

  first_part <- 0
  second_part <- 0
  two_log_two <- 2 * log(2)
  if (is.null(sigma)) { # Case where sigma is unknown
    idxs <- data %>%
      purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(8 * k0 - 6)]))
    a <- data %>%
      purrr::map2_dbl(idxs, ~ theta(.x$x, k = 4 * k0 - 3, idx = .y)) %>%
      mean()
    b <- data %>%
      purrr::map2_dbl(idxs, ~ theta(.x$x, k = 2 * k0 - 1, idx = .y)) %>%
      mean()
    c <- data %>%
      purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
      mean()
    if ((a - b > 0) & (b - c > 0) & (a - 2 * b + c > 0)) {
      first_part <- log(a - b)
      second_part <- log(b - c)
    }
  } else { # Case where sigma is known
    idxs <- data %>%
      purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(4 * k0 - 2)]))
    a <- data %>%
      purrr::map2_dbl(idxs, ~ theta(.x$x, k = 2 * k0 - 1, idx = .y)) %>%
      mean()
    b <- data %>%
      purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
      mean()
    if ((a - 2 * sigma**2 > 0) & (b - 2 * sigma**2 > 0) & (a - b > 0)) {
      first_part <- log(a - 2 * sigma**2)
      second_part <- log(b - 2 * sigma**2)
    }
  }

  (first_part - second_part) / two_log_two
}

#' Perform the estimation of H0 over a list of t0.
#'
#' @importFrom magrittr %>%
#'
#' @param data List of curves to estimate by kernel regression.
#' @param t0_list Starting times for the estimation of H0. We consider the 8k0-7
#'  nearest points of t0 for the estimation of H0 when sigma is unknown.
#' @param k0_list For the computation of the gap between the different observations.
#' @param sigma True value of sigma.
#'  If null, change estimate.
#'
#' @return A list containing the estimation of H0 at each t0.
#' @export
estimate_H0_list <- function(data, t0_list, k0_list = 2, sigma = NULL) {
  t0_list %>%
    purrr::map2_dbl(k0_list, ~ estimate_H0(data, t0 = .x, k0 = .y, sigma = sigma))
}