################################################################################
#                       Generate fractional Brownian motion                    #
################################################################################

# Load packages
library(tidyverse)

# Define functions

#' Generate fractional Brownian motion with a random noise.
#' 
#' @param M Expected number of points in the trajectory.
#'   The number of points follows a Poisson distribution with mean M.
#' @param sigma Standard deviation of the noise to add to the trajectory, can 
#'   be a list.
#' @param H Hurst coefficient.
#' 
#' @return A tibble containing the trajectory and the sampling points.
fractional_brownian_trajectory <- function(M, H, sigma, int){
  require(somebm); require(pracma)
  
  M_n <- rpois(1, M)
  t <- seq(0, 1, length.out = M_n + 1)
  x <- as.vector(fbm(hurst = H, n = M_n))
  
  if (int == TRUE) {
    x <- cumtrapz(t, x)[, 1]
    x <- cumtrapz(t, x)[, 1]
  }
  
  # Start to fill the data
  simu <- matrix(rep(0, (M_n + 1) * (length(sigma) + 2)), nrow = M_n + 1)
  simu[, 1] <- t
  simu[, 2] <- x
  
  e <- rnorm(M_n + 1, mean = 0, sd = 1)
  
  # Add columns with homoscedastic noise.
  j <- 3
  for (i in sigma) {
    simu[, j] <- x + i * e
    j = j + 1
  }

  return(as_tibble(simu, .name_repair = 'unique'))
}


#' Compute the true mean of the generated fractional Brownian motion.
#' 
#' @param t A vector of sampling points
#' @return A tibble containing the true mean trajectory and the sampling points. 
true_mean <- function(t){
  x <- rep(0, length(t))
  return(tibble(t = t, x = x))
}

#' Compute the true covariance of the generated fractional Brownian motion.
#' 
#' @param s_ A vector of sampling points
#' @param t_ A vector of sampling points
#' @param H Hurst coefficient
#' @return A tibble containing the true covariance for each pair (s, t)
true_covariance <- function(s_, t_, H){
  res <- tibble(
    s = rep(s_, times = length(t_)),
    t = rep(t_, each = length(s_))) %>% 
    mutate(phi = (abs(s)**(2*H) + abs(t)**(2*H) - abs(t - s)*(2*H))/2)
  return(res)
}

# Define some parameters
N <- 500000 # Number of curves
M <- c(1000)  # Number of points per curves (do it with 50, 200, 1000)
H <- c(0.5) # Hurst coefficient (do it with 0.4, 0.5, 0.6)
sigma <- c(0.01, 0.05, 0.1, 0.25) # Standard deviation of the noise
int <- TRUE

# Do simulation
for (m in 1:length(M)) {
    t <- seq(0, 1, length.out = M[m] + 1) # Design points
    
    simulation_ <- rerun(N, fractional_brownian_trajectory(M[m], H, sigma, int))
    mean_ <- true_mean(t)
    covariance_ <- true_covariance(t, t, H)
    
    fractional.brownian.trajectories <- list(
      simulation = simulation_,
      mean = mean_,
      covariance = covariance_,
      sigma = sigma
    )
    
    # Naming convention (fraction.brownian.trajectories-M-H-sigma)
    saveRDS(fractional.brownian.trajectories, 
            file = paste0('./data/fractional.brownian.trajectories-', 
                          M[m], '-', H, '-', paste(sigma, collapse = '.'), '.rds'))
}

