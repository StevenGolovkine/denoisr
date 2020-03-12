################################################################################
#                       Generate fractional Brownian motion                    #
################################################################################

#' Generate fractional Brownian motion with a random noise
#'
#' This function generates a realization of a fractional Brownian motion with 
#' random noise. The increments of a fractional Brownian motion are not 
#' indepenent. A fractional Brownian motion is characterized by a parameter 
#' \eqn{H}, named Hurst coefficient. We define it on \eqn{[0, 1]}.
#' 
#' @param M An integer, expected number of points in the trajectory.
#'   The number of points follows a Poisson distribution with mean \eqn{M}.
#' @param H Numeric, Hurst coefficient. \eqn{0 < H < 1}
#' @param sigma A vector of numerics, standard deviation of the noise to add to
#'  the fractional Brownian motion. 
#' @param pdf Function, probability density function for the sampling points.
#' @param L Numeric, multiplicative constant.
#'
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1} The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard deviation \eqn{\sigma}
#'  }
#'  
#' @examples 
#' fractional_brownian_trajectory(100, 0.7, 0.1)
#' fractional_brownian_trajectory(100, 0.4, 0.05, pdf = rnorm)
fractional_brownian_trajectory <- function(M, H, sigma, pdf = NULL, L = 1) {
  M_n <- rpois(1, M)
  if (!inherits(pdf, "function")) {
    t <- seq(0, 1, length.out = M_n + 1)
  } else {
    t <- pdf(M_n + 1)
    t <- t[order(t)]
  }
  x <- L * as.vector(somebm::fbm(hurst = H, n = M_n))

  # Start to fill the data
  simu <- matrix(rep(0, (M_n + 1) * (length(sigma) + 2)), nrow = M_n + 1)
  simu[, 1] <- t
  simu[, 2] <- x

  e <- rnorm(M_n + 1, mean = 0, sd = 1)

  # Add columns with homoscedastic noise.
  j <- 3
  for (i in sigma) {
    simu[, j] <- x + i * e
    j <- j + 1
  }

  dplyr::as_tibble(simu, .name_repair = "unique")
}


# Define some parameters
N <- 1000 # Number of curves
M <- 300 # Number of points per curves
H <- 0.5 # Hurst coefficient
sigma <- c(0.05) # Standard deviation of the noise


# Do simulation
t <- seq(0, 1, length.out = M + 1) # Design points

simulation_ <- purrr::rerun(N, fractional_brownian_trajectory(M, H, sigma))

fractional_brownian <- simulation_ %>% map(~ list(t = .x$...1, 
                                                  x = .x$...3, 
                                                  x_true = .x$...2))
usethis::use_data(fractional_brownian, overwrite = TRUE)
