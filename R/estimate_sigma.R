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
#'  
#' @return Numeric, an estimation of the standard deviation of the noise.
#' @export
#' @examples 
#' estimate_sigma(denoisr::fractional_brownian)
estimate_sigma <- function(data) {
  if(!inherits(data, 'list')){
    data <- checkData(data)
  }
  
  estimateSigma(data)
}
