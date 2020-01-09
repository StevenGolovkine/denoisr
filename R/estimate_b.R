################################################################################
#        Functions for bandwith parameter estimation using regularity          #
################################################################################


#' Perform an estimation of the bandwith
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
#'  \item epanechnikov (default)
#'  \item beta
#'  \item uniform
#'  }
#'
#' @return Numeric, an estimation of the bandwidth.
estimate_b <- function(data, H0 = 0.5, L0 = 1, sigma = 0, K = "epanechnikov") {

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

  # Estimate mu
  mu_hat <- data %>%
    purrr::map_int(~ length(.x$t)) %>%
    mean()

  # Estimate b
  nume <- sigma**2 * K_norm2 * factorial(floor(H0))
  deno <- H0 * L0 * K_phi
  frac <- nume / deno

  (frac / mu_hat)**(1 / (2 * H0 + 1))
}

#' Perform an estimation of the bandwidth given a list of \eqn{H_0} and \eqn{L_0}
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
#' @param sigma Numeric, an estimation of \eqn{\sigma}.
#' @param K Character string, the kernel used for the estimation:
#'  \itemize{
#'  \item epanechnikov (default)
#'  \item beta
#'  \item uniform
#'  }
#'
#' @return A vector of numerices, estimations of the bandwidth.
#' @export
#' @examples 
#' estimate_b_list(SmoothCurves::fractional_brownian, 
#'                 H0_list = 0.5, L0_list = 1, sigma = 0.05)
#' estimate_b_list(SmoothCurves::piecewise_fractional_brownian,
#'                 H0_list = c(0.2, 0.5, 0.8), L0_list = c(1, 1, 1),
#'                 sigma = 0.1, K = 'epanechnikov')
estimate_b_list <- function(data, H0_list, L0_list,
                            sigma = 0, K = "epanechnikov") {
  if (length(H0_list) != length(L0_list)) {
    stop("H0_list and L0_list must have the same length.")
  }

  H0_list %>% purrr::map2_dbl(L0_list, ~ estimate_b(data,
    H0 = .x, L0 = .y,
    sigma = sigma, K = K
  ))
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
#' estimate_b_cv(SmoothCurves::fractional_brownian)
estimate_b_cv <- function(data) {
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

  mean(unlist(bw_list))
}
