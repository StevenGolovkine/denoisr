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
#'  the piecewise fractional Brownian motion.
#'  
#' @return A tibble containing the following elements:
#'  \itemize{
#'   \item \strong{...1}: The sampling points
#'   \item \strong{...2} The true trajectory
#'   \item \strong{...3} The trajectory contaminated by noise with standard deviation \eqn{\sigma}
#'  }
#'  
#' @examples 
#' piecewise_fractional_brownian_trajectory(100, c(0.2, 0.5, 0.8), 0.1)
piecewise_fractional_brownian_trajectory <- function(M, H, sigma){
  
  M_n <- rpois(1, M)
  
  M_nn <- vector(length = length(H))
  for(i in 1:(length(H) - 1)){
    M_nn[i] <- round(M_n / length(H))
  }
  M_nn[length(H)] <- M_n - sum(M_nn)
  
  t <- c()
  for(i in 1:length(H)){
    t <- c(t, seq((i - 1) / length(H), i / length(H), length.out = M_nn[i]))
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
  simu <- matrix(rep(0, length(y) * (length(sigma) + 2)), nrow = length(y))
  simu[, 1] <- t
  simu[, 2] <- y
  
  e <- rnorm(length(y), mean = 0, sd = 1)
  
  # Add columns with homoscedastic noise.
  j <- 3
  for (i in sigma) {
    simu[, j] <- y + i * e
    j = j + 1
  }
  
  dplyr::as_tibble(simu, .name_repair = 'unique')
}


# Define some parameters
N <- 1000 # Number of curves
M <- 200  # Number of points per curves
H <- c(0.2, 0.5, 0.8) # Hurst coefficient 
sigma <- c(0.01, 0.05, 0.1) # Standard deviation of the noise 

# Do simulation
simulation_ <- purrr::rerun(N, piecewise_fractional_brownian_trajectory(M[m], H, sigma))
      
piecewise_fractional_brownian <- list(
  simulation = simulation_
)

usethis::use_data(piecewise_fractional_brownian)

