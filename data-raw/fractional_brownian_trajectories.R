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


#' Compute the mean of a fractional Brownian motion
#'
#' This function computes the mean of a fractional Brownian motion. Assume that 
#' the fractional Brownian motion is defined on \eqn{[0, 1]}, then it has expectation
#' zero for all \eqn{t} in \eqn{[0, 1]}.
#' 
#' @param t A vector of numerics, sampling points.
#' 
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item t: The sampling points
#'   \item x: The true mean trajectory
#'  }
#'  
#' @examples 
#' true_mean(seq(0, 1, length.out = 10))
true_mean <- function(t) {
  dplyr::tibble(t = t, x = rep(0, length(t)))
}

#' Compute the covariance of a fractional Brownian motion
#' 
#' This function computes the covariance of a fractional Brownian motion. Assume
#' that the fractional brownian of defined on \eqn{[0, 1]}, then it has the following
#' covariance function:
#' \deqn{E[B(t)B(s)] = \frac{1}{2}\left(\|t\|^{2H} + \|s\|^{2H} - \|t - s\|^{2H} +\right).}
#'
#' @importFrom magrittr %>%
#'
#' @param s_ A vector of numerics, sampling points.
#' @param t_ A vector of numerics, sampling points.
#' @param H Numeric, Hurst coefficient
#' 
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item s: The sampling points
#'   \item t: The sampling points
#'   \item phi: The covariance at time (s, t)
#'  }
#'  
#' @examples 
#' true_covariance(seq(0, 1, length.out = 10), seq(0, 1, length.out = 10), 0.7)
true_covariance <- function(s_, t_, H) {
  dplyr::tibble(s = rep(s_, times = length(t_)), t = rep(t_, each = length(s_))) %>%
    dplyr::mutate(phi = (abs(s)**(2 * H) + abs(t)**(2 * H) - abs(t - s) * (2 * H)) / 2)
}

# Define some parameters
N <- 1000 # Number of curves
M <- 200 # Number of points per curves
H <- 0.5 # Hurst coefficient
sigma <- c(0.01, 0.05, 0.1) # Standard deviation of the noise


# Do simulation
t <- seq(0, 1, length.out = M[m] + 1) # Design points

simulation_ <- purrr::rerun(N, fractional_brownian_trajectory(M, H, sigma))
mean_ <- true_mean(t)
covariance_ <- true_covariance(t, t, H)

fractional_brownian <- list(
  simulation = simulation_,
  mean = mean_,
  covariance = covariance_,
  sigma = sigma
)

usethis::use_data(fractional_brownian)
