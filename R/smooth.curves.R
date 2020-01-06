################################################################################
#        Functions that performs kernel smoothing over a set of curves         #
################################################################################

#' Perform a non-parametric smoothing of a set of curves.
#'
#' @importFrom magrittr %>%
#' 
#' @param data A list of curves to smooth. Each entry of the list should have
#'  two elements:
#'   - $t which correspond to the time we observed the curve.
#'   - $x which correspond to the values observed.
#' @param U Vector of points for the estimation.
#'  If NULL, smooth the curves on the same grid than they are observed.
#' @param t0 The starting time for the estimation of H0. We consider the 8k0 - 7
#'  nearest points of t0 for the estimation of H0 when sigma is unknown. Can be
#'  a list.
#' @param k0 For the computation of the gap between the different observations.
#'  Should be set as k0 = M / max(8, log(M)).
#' @param K Kernel used for the estimation.
#'  - epanechnikov (default)
#'
#' @return A list of same size of `data` containing the smoothed curves.
#' @export
smooth.curves <- function(data, U = NULL, t0 = 0.5, k0 = 2, K = "epanechnikov") {

  # Estimation of the noise
  sigma_estim <- estimate.sigma(data)

  # Estimation of H0
  H0_estim <- estimate.H0.list(data, t0_list = t0, k0_list = k0, sigma = NULL)

  # Estimation of L0
  L0_estim <- estimate.L0.list(data,
    t0_list = t0, H0_list = H0_estim,
    k0 = k0, sigma = NULL, density = NULL
  )

  # Estimation of the bandwidth
  b_estim <- estimate.b.list(data,
    H0_list = H0_estim, L0_list = L0_estim,
    sigma = sigma_estim, K = K
  )

  # Estimation of the curves
  if (is.null(U)) {
    curves <- data %>% purrr::map(~ estimate.curve(.x,
      U = .x$t, b = b_estim,
      t0_list = t0, kernel = K
    ))
  } else {
    curves <- data %>% purrr::map(~ estimate.curve(.x,
      U = U, b = b_estim,
      t0_list = t0, kernel = K
    ))
  }

  return(list(
    "parameter" = list(
      "sigma" = sigma_estim,
      "H0" = H0_estim,
      "L0" = L0_estim,
      "b" = b_estim
    ),
    "smooth" = curves
  ))
}


#' Perform a non-parametric smoothing of a set of curves when the regularity is
#' larger than 1.
#'
#' @importFrom magrittr %>%
#' 
#' @param data A list of curves to smooth. Each entry of the list should have
#'  two elements:
#'   - $t which corresponds to the time we observed the curve.
#'   - $x which corresponds to the values observed.
#'  We assume that the underlying regularity of the curve is larger than 1.
#' @param U Vector of points for the estimation.
#'  If NULL, smooth the curves on the same grid than they are observed.
#' @param t0 The starting time for the estimation of H0. We consider the 8k0 - 7
#'  nearest points of t0 for the estimation of H0 when sigma is unknown.
#' @param k0 For the computation of the gap between the different observations.
#'  Shouldd be set as k0 = M / max(8, log(M)).
#' @param K Kernel used for the estimation.
#'  - epanechnikov (default)
#' @param eps Precision parameter.
#'
#' @return A list with two entries:
#'  - parameter which contains the estimation of sigma, H0, L0 and b.
#'  - smooth which is a list of the same size than `data` containing the
#'  smoothed curves.
#' @export
smooth.curves.regularity <- function(data, U = NULL, t0 = 0.5, k0 = 2,
                                     K = "epanechnikov", eps = 0.1) {

  # Estimation of the noise
  sigma_estim <- estimate.sigma(data)

  # Estimation of H0
  H0_estim <- estimate.H0(data, t0 = t0, k0 = k0, sigma = NULL) # H > 1
  cpt <- 0
  while (H0_estim > 1 + eps) {
    L0 <- estimate.L0(data, t0 = t0, H0 = cpt + H0_estim, k0 = k0)
    b <- estimate.b(data, sigma = sigma_estim, H0 = H0_estim + cpt, L0 = L0)

    smooth <- data %>% purrr::map(~ list(
      t = .x$t,
      x = KernSmooth::locpoly(.x$t, .x$x,
        drv = 1 + cpt,
        bandwidth = b, gridsize = length(.x$t)
      )$y
    ))
    H0_estim <- estimate.H0(smooth, t0 = t0, k0 = k0, sigma = NULL)
    cpt <- cpt + 1
  }

  # Estimation of L0
  L0_estim <- estimate.L0(data, t0 = t0, H0 = cpt + H0_estim, k0 = k0)

  # Estimation of the bandwidth
  b_estim <- estimate.b(data, sigma = sigma_estim, H0 = H0_estim + cpt, L0 = L0_estim)

  # Estimation of the curves
  if (is.null(U)) {
    curves <- data %>% purrr::map(~ estimate.curve(.x,
      U = .x$t, b = b_estim,
      t0_list = t0, kernel = K
    ))
  } else {
    curves <- data %>% purrr::map(~ estimate.curve(.x,
      U = U, b = b_estim,
      t0_list = t0, kernel = K
    ))
  }

  return(list(
    "parameter" = list(
      "sigma" = sigma_estim,
      "H0" = H0_estim + cpt,
      "L0" = L0_estim,
      "b" = b_estim
    ),
    "smooth" = curves
  ))
}
