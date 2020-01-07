######################################################################################
#           Functions for bandwith parameter estimation using regularity             #
######################################################################################


#' Perform the estimation of the bandwith
#'
#' @importFrom magrittr %>%
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param H0 An estimation of H0.
#' @param L0 An estimation of L0.
#' @param sigma An estimation of sigma.
#' @param K Which kernel to use?
#'  - epanechnikov (default)
#'  - beta
#'  - uniform
#'
#' @return An estimation of H0.
estimate.b <- function(data, H0 = 0.5, L0 = 1, sigma = 0, K = "epanechnikov") {
  S_N <- data

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
  mu_hat <- S_N %>%
    purrr::map_int(~ length(.x$t)) %>%
    mean()

  # Estimate b
  # H0_tilde <- H0 - (1 / log(mu_hat))
  nume <- (sigma**2 * K_norm2 * factorial(floor(H0)))
  deno <- (H0 * L0 * K_phi)
  frac <- nume / deno
  b_hat <- (frac / mu_hat)**(1 / (2 * H0 + 1))

  return(b_hat)
}

#' Perform the estimation of the bandwidth over a list of H0 and t0.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' 
#' @param data List of curves to estimate by kernel regression.
#' @param H0_list List of estimates of H0.
#' @param L0_list <- List of estimates of L0.
#' @param sigma An estimation of sigma.
#' @param K Which kernel to use?
#'  - epanechnikov
#'  - beta
#'  - uniform
#'
#' @return List of estimates of the bandwidth.
#' @export
estimate_b_list <- function(data, H0_list, L0_list,
                            sigma = 0, K = "epanechnikov") {
  if (length(H0_list) != length(L0_list)) {
    stop("H0_list and L0_list must have the same length.")
  }

  b_hat_list <- H0_list %>%
    purrr::map2_dbl(L0_list, ~ estimate.b(data,
      H0 = .x, L0 = .y,
      sigma = sigma, K = K
    ))
  return(b_hat_list)
}

#' Perform the estimation of the bandwith using CV.
#'
#' @importFrom magrittr %>%
#'
#' @param data List of curves to estimate by kernel regression.
#'
#' @return An estimation of the bandwidth by cross-validation.
#' @export
estimate.b.cv <- function(data) {

  cl <- parallel::detectCores() %>%
    -1 %>%
    parallel::makeCluster()
  doParallel::registerDoParallel(cl)

  j = 1:length(data)
  bw_list <- foreach(j = iterators::iter(j)) %dopar% {
    sqrt(5) * np::npregbw(x ~ t,
      data = data[[j]],
      bwmethod = "cv.ls", # Least Square Cross Validation
      ckertype = "epanechnikov", # Kernel used
      regtype = "lc" # Local Constant Regression
    )$bw 
  }

  parallel::stopCluster(cl)

  return(mean(unlist(bw_list)))
}
