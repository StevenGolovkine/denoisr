################################################################################
#        Functions that performs kernel smoothing over a set of curves         #
################################################################################

#' Perform a non-parametric smoothing of a set of curves
#'
#' This function performs a non-parametric smoothing of a set of curves using the
#' Nadaraya-Watson estimator. The bandwidth is estimated using the method from 
#' \cite{add ref}.
#' 
#' @importFrom magrittr %>%
#'
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U A vector of numerics, sampling points at which estimate the curves.
#'  If NULL, the sampling points for the estimation are the same than the 
#'  observed ones.
#' @param t0_list A vector of numerics, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0_list A vector of numerics, the number of neighbors of \eqn{t_0} to 
#'  consider. Should be set as \deqn{k0 = (M / log(M) + 7) / 8}. We can set a 
#'  different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just 
#'  put a unique numeric.
#' @param K Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item uniform
#'   \item beta
#'  }
#'
#' @return A list, which contains two elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} An estimation of the standard deviation of the noise
#'   \item \strong{H0} An estimation of \eqn{H_0}
#'   \item \strong{L0} An estimation of \eqn{L_0}
#'   \item \strong{b} An estimation of the bandwidth
#'  }
#'  The second one is another list which contains the estimation of the curves:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  } 
#' @export
#' @examples 
#' df <- smooth_curves(SmoothCurves::fractional_brownian)
#' df <- smooth_curves(SmoothCurves::piecewise_fractional_brownian, 
#'                     t0_list = c(0.15, 0.5, 0.85), k0_list = 6)
smooth_curves <- function(data, U = NULL, 
                          t0_list = 0.5, k0_list = 2, K = "epanechnikov") {

  if(!inherits(data, 'list')){
    data <- checkData(data)
  }
  
  # Estimation of the noise
  sigma_estim <- estimate_sigma(data)

  # Estimation of H0
  H0_estim <- estimate_H0_list(data, t0_list = t0_list, k0_list = k0_list, sigma = NULL)

  # Estimation of L0
  L0_estim <- estimate_L0_list(data,
    t0_list = t0_list, H0_list = H0_estim,
    k0 = k0_list[1], sigma = NULL, density = FALSE
  )

  # Estimation of the bandwidth
  b_estim <- estimate_b_list(data,
    H0_list = H0_estim, L0_list = L0_estim,
    sigma = sigma_estim, K = K
  )

  # Estimation of the curves
  if (is.null(U)) {
    curves <- data %>% purrr::map(~ estimate_curve(.x,
      U = .x$t, b = b_estim,
      t0_list = t0_list, kernel = K
    ))
  } else {
    curves <- data %>% purrr::map(~ estimate_curve(.x,
      U = U, b = b_estim,
      t0_list = t0_list, kernel = K
    ))
  }

  list(
    "parameter" = list(
      "sigma" = sigma_estim,
      "H0" = H0_estim,
      "L0" = L0_estim,
      "b" = b_estim
    ),
    "smooth" = curves
  )
}


#' Perform a non-parametric smoothing of a set of curves when the regularity is
#' larger than 1
#'
#' This function performs a non-parametric smoothing of a set of curves using the
#' Nadaraya-Watson estimator when the regularity of the underlying curves is
#' larger than 1. The bandwidth is estimated using the method from 
#' \cite{add ref}. In the case of a regularly larger than 1, we currently 
#' assume that the regularly is the same all over the curve.
#' 
#' @importFrom magrittr %>%
#'
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U A vector of numerics, sampling points at which estimate the curves.
#'  If NULL, the sampling points for the estimation are the same than the 
#'  observed ones.
#' @param t0 Numeric, the sampling point at which we estimate \eqn{H0}. We will 
#'  consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for the estimation of
#'  \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0 Numeric, the number of neighbors of \eqn{t_0} to consider. Should 
#'  be set as \eqn{k0 = (M / log(M) + 7) / 8}.
#' @param K Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item uniform
#'   \item beta
#'  }
#' @param eps Numeric, precision parameter. It is used to control how much larger 
#'  than 1, we have to be in order to consider to have a regularity larger than 1
#'  (default to 0.1).
#'
#' @return A list, which contains two elements. The first one is a list which 
#'  contains the estimated parameters:
#'  \itemize{
#'   \item \strong{sigma} An estimation of the standard deviation of the noise
#'   \item \strong{H0} An estimation of \eqn{H_0}
#'   \item \strong{L0} An estimation of \eqn{L_0}
#'   \item \strong{b} An estimation of the bandwidth
#'  }
#'  The second one is another list which contains the estimation of the curves:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The estimated points.
#'  } 
#' @export
#' @examples 
#' df <- smooth_curves_regularity(SmoothCurves::fractional_brownian)
smooth_curves_regularity <- function(data, U = NULL, t0 = 0.5, k0 = 2,
                                     K = "epanechnikov", eps = 0.1) {

  if(!inherits(data, 'list')){
    data <- checkData(data)
  }
  
  # Estimation of the noise
  sigma_estim <- estimate_sigma(data)

  # Estimation of H0
  H0_estim <- estimate_H0(data, t0 = t0, k0 = k0, sigma = NULL) # H > 1
  cpt <- 0
  while (H0_estim > 1 + eps) {
    L0 <- estimate_L0(data, t0 = t0, H0 = cpt + H0_estim, k0 = k0)
    b <- estimate_b(data, sigma = sigma_estim, H0 = H0_estim + cpt, L0 = L0)

    smooth <- data %>% purrr::map(~ list(
      t = .x$t,
      x = KernSmooth::locpoly(.x$t, .x$x,
        drv = 1 + cpt,
        bandwidth = b, gridsize = length(.x$t)
      )$y
    ))
    H0_estim <- estimate_H0(smooth, t0 = t0, k0 = k0, sigma = NULL)
    cpt <- cpt + 1
  }

  # Estimation of L0
  L0_estim <- estimate_L0(data, t0 = t0, H0 = cpt + H0_estim, k0 = k0)

  # Estimation of the bandwidth
  b_estim <- estimate_b(data, sigma = sigma_estim, H0 = H0_estim + cpt, L0 = L0_estim)

  # Estimation of the curves
  if (is.null(U)) {
    curves <- data %>% purrr::map(~ estimate_curve(.x,
      U = .x$t, b = b_estim,
      t0_list = t0, kernel = K
    ))
  } else {
    curves <- data %>% purrr::map(~ estimate_curve(.x,
      U = U, b = b_estim,
      t0_list = t0, kernel = K
    ))
  }

  list(
    "parameter" = list(
      "sigma" = sigma_estim,
      "H0" = H0_estim + cpt,
      "L0" = L0_estim,
      "b" = b_estim
    ),
    "smooth" = curves
  )
}
