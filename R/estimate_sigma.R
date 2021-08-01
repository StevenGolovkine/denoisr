################################################################################
#       Functions for sigma parameter estimation using regularity              #
################################################################################


#' Perform an estimation of the standard deviation of the noise
#' 
#' This function performs an estimation of the standard deviation of the noise 
#' in the curves. Based on \cite{add ref}, the following formula is used:
#' \deqn{\hat{\sigma^2} = \frac{1}{N}\sum_{n = 1}^{N} 
#'       \frac{1}{2(M_n - 1)}\sum_{l = 2}^{M_n}(Y_{n, (l)} - Y_{n, (l-1)})^2}
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param t0 Numeric, the sampling point at which we estimate \eqn{\sigma}. We will 
#'  consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for the estimation of
#'  \eqn{\sigma}
#' @param k0 Numeric, the number of neighbors of \eqn{t_0} to consider. Should 
#'  be set as \eqn{k0 = M * exp(-log(log(M))^2)}. 
#'  
#' @return Numeric, an estimation of the standard deviation of the noise.
#' @export
#' @examples
#' df <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
#' estimate_sigma(df, t0 = 0.5, k0 = 14)
estimate_sigma <- function(data, t0 = 0.5, k0 = 2){
  if(!inherits(data, 'list')) data <- checkData(data)
  idxs <- data %>% purrr::map_dbl(~ min(order(abs(.x$t - t0))[seq_len(8 * k0 - 6)]))
  df_sub <- data %>% purrr::map2(idxs, ~ list(t = .x$t[.y:(.y + 8 * k0 - 7)],
                                              x = .x$x[.y:(.y + 8 * k0 - 7)]))
  estimateSigma(df_sub)
}


#' Perform an estimation of the noise given a list of \eqn{t_0}.
#' 
#' This function performs an estimation of the standard deviation of the noise
#' at different times. 
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @param t0_list A vector of numeric, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{H_0} when \eqn{\sigma} is unknown.
#' @param k0_list A vector of numeric, the number of neighbors of \eqn{t_0} to 
#'  consider. Should be set as \deqn{k0 = M * exp(-(log(log(M))**2))}. We can set a 
#'  different \eqn{k_0}, but in order to use the same for each \eqn{t_0}, just 
#'  put a unique numeric.
#'  
#' @return A list, an estimation of sigma at different \eqn{t_0}.
#' @export
#' @examples
#' df <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
#' estimate_sigma_list(df, t0_list = c(0.25, 0.5, 0.75), k0 = 14)
estimate_sigma_list <- function(data, t0_list, k0_list){
  if(!inherits(data, 'list')) data <- checkData(data)
  t0_list %>%
    purrr::map2_dbl(k0_list, ~ estimate_sigma(data, t0 = .x, k0 = .y))
}
