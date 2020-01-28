###############################################################################
#        Functions for L0 parameter estimation using regularity               #
###############################################################################


#' Perform an estimation of \eqn{L_0}
#'
#' This function performs an estimation of \eqn{H_0} used for the estimation of
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data using the method of \cite{add ref}. 
#' 
#' @importFrom magrittr %>%
#'
#' @family estimate \eqn{L_0}
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param t0 Numeric, the sampling point at which we estimate \eqn{H0}. We will 
#'  consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for the estimation of
#'  \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param H0 Numeric, an estimation of \eqn{H_0}
#' @param k0 Numeric, the number of neighbors of \eqn{t_0} to consider. Should 
#'  be set as \eqn{k0 = (M / log(M) + 7) / 8}.
#' @param sigma Numeric, true value of sigma. Can be NULL.
#' @param density Logical, do the sampling points have a uniform distribution? 
#'  (default is FALSE)
#'
#' @return Numeric, an estimation of L0.
estimate_L0 <- function(data, t0 = 0, H0 = 0,
                        k0 = 2, sigma = NULL, density = FALSE) {

  # Estimate mu
  mu_hat <- data %>%
    purrr::map_int(~ length(.x$t)) %>%
    mean()

  # Estimate L0
  theta <- function(v, k, idx) (v[idx + 2 * k - 1] - v[idx + k])**2
  eta <- function(v, k, idx, H) (v[idx + 2 * k - 1] - v[idx + k])**(2 * H)

  nume <- 1
  deno <- 1
  if (!density) { # Case where the density is not known
    if (is.null(sigma)) { # Subcase where sigma is not known
      idxs <- data %>%
        purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(4 * k0 - 2)]))
      a <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = 2 * k0 - 1, idx = .y)) %>%
        mean()
      b <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
        mean()
      c <- data %>%
        purrr::map2_dbl(idxs, ~ eta(.x$t, k = 2 * k0 - 1, idx = .y, H = H0)) %>%
        mean()
      d <- data %>%
        purrr::map2_dbl(idxs, ~ eta(.x$t, k = k0, idx = .y, H = H0)) %>%
        mean()
      if ((a - b > 0) & (c - d > 0)) {
        nume <- a - b
        deno <- c - d
      }
    } else { # Subcase where sigma is known
      idxs <- data %>%
        purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(2 * k0)]))
      a <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
        mean()
      b <- data %>%
        purrr::map2_dbl(idxs, ~ eta(.x$t, k = k0, idx = .y, H = H0)) %>%
        mean()
      if ((a - 2 * sigma**2 > 0) & b > 0) {
        nume <- a - 2 * sigma**2
        deno <- b
      }
    }
  } else { # Case where the density is known (only the uniform case)
    if (is.null(sigma)) { # Subcase where sigma is not known
      idxs <- data %>%
        purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(4 * k0 - 2)]))
      a <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = 2 * k0 - 1, idx = .y)) %>%
        mean()
      b <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
        mean()
      if (a - b > 0) nume <- a - b
      deno <- (2**(2 * H0) - 1) * ((k0 - 1) / (mu_hat + 1))**(2 * H0)
    } else { # Subcase where sigma is known
      idxs <- data %>%
        purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(2 * k0)]))
      a <- data %>%
        purrr::map2_dbl(idxs, ~ theta(.x$x, k = k0, idx = .y)) %>%
        mean()
      if (a - 2 * sigma**2 > 0) nume <- a - 2 * sigma**2
      deno <- ((k0 - 1) / (mu_hat + 1))**(2 * H0)
    }
  }

  (nume / deno)**0.5
}


#' Perform the estimation of \eqn{L_0} given a list of \eqn{t_0}
#'
#' This function performs an estimation of \eqn{H_0} used for the estimation of
#' the bandwidth for a univariate kernel regression estimator defined over 
#' continuous domains data using the method of \cite{add ref}. 
#'
#' @importFrom magrittr %>%
#' 
#' @family estimate \eqn{L_0}
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param t0_list A vetor of numerics, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param H0_list A vector of numerics, an estimation of \eqn{H_0} at every 
#'  \eqn{t_0} given in \code{t0_list}.
#' @param k0 Numeric, the number of neighbors of \eqn{t_0} to consider. Should 
#'  be set as \eqn{k0 = (M / log(M) + 7) / 8}.
#' @param sigma Numeric, true value of sigma. Can be NULL.
#' @param density Logical, do the sampling points have a uniform distribution? 
#'  (default is FALSE)
#'
#' @return A vector of numerics, an estimation of \eqn{L_0} at each \eqn{t_0}.
#' @export
#' @examples 
#' estimate_L0_list(SmoothCurves::fractional_brownian, 
#'                 t0_list = 0.5, H0_list = 0.5)
#' estimate_L0_list(SmoothCurves::piecewise_fractional_brownian,
#'                 t0_list = c(0.15, 0.5, 0.85), H0_list = c(0.2, 0.5, 0.8), 
#'                 k0 = 6)
estimate_L0_list <- function(data, t0_list, H0_list,
                             k0 = 2, sigma = NULL, density = FALSE) {
  if(!inherits(data, 'list')){
    data <- checkData(data)
  }
  
  if (length(t0_list) != length(H0_list)) {
    stop("t0_list and H0_list must have the same length")
  }

  t0_list %>%
    purrr::map2_dbl(H0_list, ~ estimate_L0(data,
      t0 = .x, H0 = .y,
      k0 = k0, sigma = sigma, density = density
    ))
}
