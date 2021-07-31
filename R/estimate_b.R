################################################################################
#        Functions for bandwidth parameter estimation using regularity         #
################################################################################

#' Perform an estimation of the bandwidth for the smoothing of curves.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data using the method of 
#' \cite{add ref}. An estimation of \eqn{H_0}, \eqn{L_0} and \eqn{\sigma} have 
#' to be provided to estimate the bandwidth. 
#'
#' @importFrom magrittr %>%
#'
#' @family estimate bandwidth
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param H0 Numeric, an estimation of \eqn{H_0}.
#' @param L0 Numeric, an estimation of \eqn{L_0}.
#' @param sigma Numeric, an estimation of \eqn{\sigma}.
#' @param K Character string, the kernel used for the estimation:
#'  \itemize{
#'   \item epanechnikov (default)
#'   \item beta
#'   \item uniform
#'  }
#'
#' @return Numeric, an estimation of the bandwidth.
#' @export
#' @examples
#' X <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
#' estimate_b(X, H0 = 0.5, L0 = 1, sigma = 0.05)
estimate_b <- function(data, H0 = 0.5, L0 = 1, sigma = 0, K = "epanechnikov"){
  if(!inherits(data, 'list')) data <- checkData(data)
  
  # Set kernel constants
  if (K == "epanechnikov") {
    K_norm2 <- 0.6
    K_phi <- 3 / ((H0 + 1) * (H0 + 3))
  } else if (K == "beta") {
    K_norm2_f <- function(x, alpha = 1, beta = 1) {
      x**(2 * (alpha - 1)) * (1 - x)**(2 * (beta - 1)) / beta(alpha, beta)**2
    }
    phi <- function(x, H0, alpha = 1, beta = 1) {
      abs(x**(alpha - 1) * (1 - x)**(beta - 1)) * x**H0 / abs(beta(alpha, beta))
    }
    K_norm2 <- stats::integrate(K_norm2_f, lower = 0, upper = 1)$value
    K_phi <- stats::integrate(phi, lower = 0, upper = 1, H0 = H0)$value
  } else {
    K_norm2 <- 1
    K_phi <- 1 / (H0 + 1)
  }

  M_n <- data %>%
    purrr::map_int(~ length(.x$t))

  nume <- sigma**2 * factorial(floor(H0))**2 * K_norm2
  deno <- 2 * L0**2 * H0 * K_phi
  frac <- nume / deno
  (frac / M_n)**(1 / (2 * H0 + 1))
}


#' Perform an estimation of the bandwidth given a list of \eqn{H_0} and \eqn{L_0}
#' for the smoothing of curves.
#' 
#' This function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data using the method of 
#' \cite{add ref}. An estimation of \eqn{H_0}, \eqn{L_0} and \eqn{\sigma} have 
#' to be provided to estimate the bandwidth. 
#'
#' @importFrom magrittr %>%
#'
#' @family estimate bandwidth
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param H0_list A vector of numerics, estimations of \eqn{H_0}.
#' @param L0_list A vector of numerics, estimations of \eqn{L_0}.
#' @param sigma A vector of numerics, an estimation of \eqn{\sigma}.
#' @param K Character string, the kernel used for the estimation:
#'  \itemize{
#'  \item epanechnikov (default)
#'  \item beta
#'  \item uniform
#'  }
#'
#' @return A vector of numerics, estimations of the bandwidth.
#' @export
#' @examples
#' X <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
#' estimate_b_list(X, H0_list = 0.5, L0_list = 1, sigma = 0.05)
#' 
#' X <- generate_piecewise_fractional_brownian(N = 1000, M = 300, 
#'                                             H = c(0.2, 0.5, 0.8), 
#'                                             sigma = 0.05)
#' estimate_b_list(X, H0_list = c(0.2, 0.5, 0.8), L0_list = c(1, 1, 1), 
#'                 sigma = 0.1, K = 'epanechnikov')
estimate_b_list <- function(data, H0_list, L0_list, sigma = 0, 
                            K = "epanechnikov") {
  if(!inherits(data, 'list')) data <- checkData(data)
  
  if (length(H0_list) != length(L0_list)) {
    stop("H0_list and L0_list must have the same length.")
  }

  purrr::pmap_dfc(list(H0_list, L0_list, sigma),
                  function(H0, L0, s){
                    estimate_b(data, H0 = H0, L0 = L0, sigma = s, K = K)
                  })
}

#' Perform an estimation of the bandwidth
#'
#' This function performs an estimation of the bandwidth to be used in the
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
#' @param t0_list A vector of numerics, the sampling points at which we estimate 
#'  \eqn{H0}. We will consider the \eqn{8k0 - 7} nearest points of \eqn{t_0} for 
#'  the estimation of \eqn{L_0} when \eqn{\sigma} is unknown.
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
#' @return A list, with elements:
#'  \itemize{
#'   \item \strong{sigma} An estimation of the standard deviation of the noise
#'   \item \strong{H0} An estimation of \eqn{H_0}
#'   \item \strong{L0} An estimation of \eqn{L_0}
#'   \item \strong{b} An estimation of the bandwidth
#'  }
#' @export
#' @examples
#' X <- generate_fractional_brownian(N = 1000, M = 300,  H = 0.5, sigma = 0.05)
#' estimate_bandwidth(X)
#' 
#' X <- generate_piecewise_fractional_brownian(N = 1000, M = 300, 
#'                                             H = c(0.2, 0.5, 0.8), 
#'                                             sigma = 0.05)
#' estimate_bandwidth(X, t0_list = c(0.15, 0.5, 0.85), k0_list = 6)
estimate_bandwidth <- function(data, t0_list = 0.5, k0_list = 2,
                               K = "epanechnikov") {
  if(!inherits(data, 'list')) data <- checkData(data)
  
  # Estimation of the noise
  sigma_estim <- estimate_sigma_list(data, t0_list, k0_list)
  
  # Estimation of H0
  H0_estim <- estimate_H0_list(data, t0_list, k0_list)

  # Estimation of L0
  L0_estim <- estimate_L0_list(data, t0_list = t0_list, H0_list = H0_estim,
                               k0_list = k0_list, sigma = NULL, density = FALSE)
  
  # Estimation of the bandwidth
  b_estim <- estimate_b_list(data, H0_list = H0_estim, L0_list = L0_estim,
                             sigma = sigma_estim, K = K)
  
  list(
    "sigma" = sigma_estim,
    "H0" = H0_estim,
    "L0" = L0_estim,
    "b" = b_estim
  )
}

#' Perform an estimation of the bandwidth using least-squares cross validation
#' 
#' The function performs an estimation of the bandwidth for a univariate kernel
#' regression estimator defined over continuous data using least-square cross 
#' validation for each curve and return the average bandwidth among them.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @family estimate bandwidth
#' @seealso \code{\link[np]{npregbw}}
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#'
#' @return Numeric, an estimation of the bandwidth.
#' @export
#' @examples 
#' X <- generate_fractional_brownian(N = 5, M = 300, H = 0.5, sigma = 0.05)
#' estimate_b_cv(X)
estimate_b_cv <- function(data) {
  if(!inherits(data, 'list')) data <- checkData(data)
  
  # Create clusters for parallel computation
  cl <- parallel::detectCores() %>%
    -1 %>%
    parallel::makeCluster()
  doParallel::registerDoParallel(cl)

  j <- 1:length(data)
  bw_list <- foreach(j = iterators::iter(j)) %dopar% {
    sqrt(5) * np::npregbw(x ~ t,
      data = data[[j]],
      bwmethod = "cv.ls", # Least Square Cross Validation
      ckertype = "epanechnikov", # Kernel used
      regtype = "lc" # Local Constant Regression
    )$bw
  }

  parallel::stopCluster(cl)

  bw_list
}
