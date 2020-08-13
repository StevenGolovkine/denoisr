################################################################################
#                   Generate fractional Brownian motion                        #
################################################################################

#' Generate fractional Brownian motion with a random noise
#'
#' This function generates a realization of a fractional Brownian motion with 
#' random noise. The increments of a fractional Brownian motion are not 
#' independent. A fractional Brownian motion is characterized by a parameter 
#' \eqn{H}, named Hurst coefficient. We define it on \eqn{[0, 1]}.
#' 
#' @importFrom stats rnorm
#' @importFrom stats rpois
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
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
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


#' Generate a list of fractional Brownian trajectory.
#'
#' This function generates a list of realizations of a fractional Brownian 
#' motion with random noise. The increments of a fractional Brownian motion 
#' are not independent. A fractional Brownian motion is characterized by a 
#' parameter \eqn{H}, named Hurst coefficient. We define it on \eqn{[0, 1]}.
#' 
#' @param N An integer, number of curves to simulate.
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
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
#'  
#' @export
#' 
#' @examples
#' generate_fractional_brownian(100, 50, 0.7, 0.1)
#' generate_fractional_brownian(100, 50, 0.5, 0.05, pdf = rnorm)
generate_fractional_brownian <- function(N = 100, M = 10, H = 0.5, 
                                         sigma = 0.05, pdf = NULL, L = 1){
  
  simulation_ <- purrr::rerun(N, fractional_brownian_trajectory(M, H, sigma, 
                                                                pdf, L))
  purrr::map(simulation_, ~ list(t = .x$...1, x = .x$...3, x_true = .x$...2))  
}


################################################################################
#               Generate piecewise fractional Brownian motion                  #
################################################################################

#' Generate piecewise fractional Brownian motion with a random noise
#' 
#' This function generates a realization of a piecewise fractional Brownian motion
#' wih random noise. A piecewise fractional Brownian motion is defined by a non
#' constant Hurst parameter along the sampling points. We observe the process at
#' regularly spaced time \eqn{t_i = \frac{i}{M_n}}, where \eqn{i = 0, \dots, M_n}.
#' We define a segmentation \eqn{\tau = (\tau_k)_{k=0, \dots, K+1}}, with 
#' \eqn{0 = \tau_0 < \tau_1 < \dots < \tau_{K} < \tau_{K+1} = 1}. So, on the 
#' interval \eqn{[\tau_k, \tau_{k+1}]}, for \eqn{k = 0, \dots, K}, the process 
#' is a fractional Brownian motion with Hurst parameter \eqn{H_k}. 
#' 
#' @param M An integer, expected number of points in the trajectory.
#'   THe number of points follows a Poisson distribution with mean \eqn{M}.
#' @param H A vector of numeric, Hurst coefficients. \eqn{0 < H_k < 1}
#' @param sigma A vector of numerics, standard deviation of the noise to add to 
#'  the piecewise fractional Brownian motion. Should have the length of H. It
#'  adds heteroscedastic noise to the data.
#' @param pdf A function for the generation of the sampling points.
#'  
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1}: The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
piecewise_fractional_brownian_trajectory <- function(M, H, sigma, pdf = NULL){
  
  M_n <- rpois(1, M)
  
  M_nn <- vector(length = length(H))
  for(i in 1:(length(H) - 1)){
    M_nn[i] <- round(M_n / length(H))
  }
  M_nn[length(H)] <- M_n - sum(M_nn)
  
  if (!inherits(pdf, "function")) {
    t <- c()
    for(i in 1:length(H)){
      t <- c(t, seq((i - 1) / length(H), i / length(H), length.out = M_nn[i]))
    }
  } else {
    t <- pdf(M_n)
    t <- t[order(t)]
  }

  x <- list()
  for(i in seq_along(H)){
    x[i] <- list(as.vector(somebm::fbm(hurst = H[i], n = M_nn[i] - 1)))
  }
  
  # Some continuity in the change point
  y <- c(x[[1]])
  for(i in 2:length(x)){
    y <- c(y, x[[i]] + y[length(y)])
  }
  
  # Start to fill the data
  simu <- matrix(rep(0, length(y) * 3), nrow = length(y))
  simu[, 1] <- t
  simu[, 2] <- y
  
  e <- rnorm(length(y), mean = 0, sd = 1)
  
  if (length(sigma) > 1){
    sigma <- rep(sigma, M_nn)
  }
  
  # Add columns with homoscedastic noise.
  simu[, 3] <- y + sigma * e
  
  dplyr::as_tibble(simu, .name_repair = 'unique')
}


#' Generate a list of piecewise fractional Brownian trajectory.
#' 
#' This function generates a list of realizations of a piecewise fractional Brownian 
#' motion with random noise. A piecewise fractional Brownian motion is defined by a 
#' non constant Hurst parameter along the sampling points. We observe the process at
#' regularly spaced time \eqn{t_i = \frac{i}{M_n}}, where \eqn{i = 0, \dots, M_n}.
#' We define a segmentation \eqn{\tau = (\tau_k)_{k=0, \dots, K+1}}, with 
#' \eqn{0 = \tau_0 < \tau_1 < \dots < \tau_{K} < \tau_{K+1} = 1}. So, on the 
#' interval \eqn{[\tau_k, \tau_{k+1}]}, for \eqn{k = 0, \dots, K}, the process 
#' is a fractional Brownian motion with Hurst parameter \eqn{H_k}. 
#' 
#' @param N An integer, number of curves to simulate.
#' @param M An integer, expected number of points in the trajectory.
#'   THe number of points follows a Poisson distribution with mean \eqn{M}.
#' @param H A vector of numeric, Hurst coefficients. \eqn{0 < H_k < 1}
#' @param sigma A vector of numerics, standard deviation of the noise to add to 
#'  the piecewise fractional Brownian motion. Should have the length of H. It
#'  adds heteroscedastic noise to the data.
#' @param pdf A function for the generation of the sampling points.
#'  
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1}: The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
#'
#' @export
#' @examples 
#' generate_piecewise_fractional_brownian(100, 50, c(0.2, 0.5, 0.8), 0.1)
generate_piecewise_fractional_brownian <- function(N = 100, M = 100, 
                                                   H = c(0.2, 0.5, 0.8), 
                                                   sigma = 0.05,
                                                   pdf = NULL){
  
  simulation_ <- purrr::rerun(N, piecewise_fractional_brownian_trajectory(M, H,
                                                                          sigma,
                                                                          pdf))
  purrr::map(simulation_, ~ list(t = .x$...1, x = .x$...3, x_true = .x$...2))
}


################################################################################
#               Generate integrate fractional Brownian motion                  #
################################################################################

#' Generate integrate fractional Brownian motion with a random noise
#'
#' This function generates a realization of an integrate fractional Brownian 
#' motion with random noise. The increments of a integrate fractional Brownian 
#' motion are not independent. An integrate fractional Brownian motion is 
#' characterized by a parameter \eqn{H}, named Hurst coefficient. 
#' We define it on \eqn{[0, 1]}.
#' 
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' 
#' @param M An integer, expected number of points in the trajectory.
#'   The number of points follows a Poisson distribution with mean \eqn{M}.
#' @param H Numeric, Hurst coefficient. \eqn{0 < H < 1}. As we return its
#'  integrated version, the true Hurst will be 1 + H.
#' @param sigma A vector of numerics, standard deviation of the noise to add to
#'  the fractional Brownian motion.
#' @param L Numeric, multiplicative constant.
#'
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1} The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
integrate_fractional_brownian_trajectory <- function(M, H, sigma, L = 1) {
  M_n <- rpois(1, M)
  t <- seq(0, 1, length.out = M_n + 1)
  
  x <- L * as.vector(somebm::fbm(hurst = H, n = M_n))
  
  # Start to fill the data
  simu <- matrix(rep(0, (M_n + 1) * (length(sigma) + 2)), nrow = M_n + 1)
  simu[, 1] <- t
  simu[, 2] <- cumsum(x) / (M_n + 1)
  
  e <- rnorm(M_n + 1, mean = 0, sd = 1)

  # Add columns with homoscedastic noise.
  j <- 3
  for (i in sigma) {
    simu[, j] <- cumsum(x) / (M_n + 1) + i * e
    j <- j + 1
  }
  
  dplyr::as_tibble(simu, .name_repair = "unique")
}


#' Generate a list of integrate fractional Brownian trajectory.
#'
#' This function generates a list of realizations of a integrate fractional 
#' Brownian motion with random noise. The increments of a integrate fractional 
#' Brownian motion are not independent. An integrate fractional Brownian motion 
#' is characterized by a parameter \eqn{H}, named Hurst coefficient. 
#' We define it on \eqn{[0, 1]}.
#' 
#' @param N An integer, number of curves to simulate.
#' @param M An integer, expected number of points in the trajectory.
#'   The number of points follows a Poisson distribution with mean \eqn{M}.
#' @param H Numeric, Hurst coefficient. \eqn{0 < H < 1}. As we return its
#'  integrated version, the true Hurst will be 1 + H.
#' @param sigma A vector of numerics, standard deviation of the noise to add to
#'  the fractional Brownian motion.
#' @param L Numeric, multiplicative constant.
#'
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1} The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard 
#'   deviation \eqn{\sigma}
#'  }
#'  
#' @export
#' @examples
#' generate_integrate_fractional_brownian(100, 50, 0.7, 0.1)
#' generate_integrate_fractional_brownian(100, 50, 0.5, 0.05, 2)
generate_integrate_fractional_brownian <- function(N = 100, M = 10, H = 0.5, 
                                                   sigma = 0.05, L = 1){
  
  simulation_ <- purrr::rerun(N, integrate_fractional_brownian_trajectory(M, H, sigma, L))
  purrr::map(simulation_, ~ list(t = .x$...1, x = .x$...3, x_true = .x$...2))  
}
