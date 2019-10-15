######################################################################################
#                          Generate quadratic trajectories                           #
######################################################################################

# Load packages
library(tidyverse)

# Define functions

#' Generate quadratic trajectory with a random noise.
#' 
#' @param M Number of points in the trajectory
#' @param sigma Standard deviation of the noise to add to the trajectory
#' @param hetero_noise Should the moise be heteroscedastic? (default = FALSE)
#' @param hetero_func A function to generate the noise
#' @return A tibble containing the trajectory and the sampling points.
quadratic_trajectory <- function(M, sigma, hetero_noise = FALSE, hetero_func = NULL){
  if(is.null(hetero_func)){
    sd.function <- function(t) (abs(t)+1)**(5/6)
  } else{
    sd.function <- hetero_func
  }
  
  t <- seq(from = 0, to = 1, length.out = M+1)
  e <- rnorm(M+1, mean = 0, sd = sigma)
  x <- (t - 0.5)**2
  if(hetero_noise) e <- sd.function((x - median(x)) / sd(x)) * e
  return(tibble(t = t, x = x + e))
}

#' Compute the true mean of the generated quadratic trajectories.
#' 
#' @param t A vector of sampling points
#' @return A tibble containing the true mean trajectory and the sampling points. 
true_mean <- function(t){
  x <- (t - 0.5)**2
  return(tibble(t = t, x = x))
}

#' Compute the true covariance of the generated quadratic trajectories.
#' 
#' @param s_ A vector of sampling points
#' @param t_ A vector of sampling points
#' @return A tibble containing the true covariance for each pair (s, t)
true_covariance <- function(s_, t_){
  phi <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
  return(tibble(
    s = rep(s_, times = length(t_)),
    t = rep(t_, each = length(s_)),
    phi = as.vector(phi)))
}


# Define some parameters
N <- 10000 # Number of curves
M <- c(50, 200, 1000)  # Number of points per curves (do it with 50, 200, 1000)
sigma <- c(0, 0.01, 0.05, 0.1) # Standard deviation of the noise (do it with 0, 0.01, 0.05, 0.1)
hetero_noise <- TRUE # Add heteroscedastic noise?
hetero_func <- function(x) return (2 * exp(abs(x)) / (1 + exp(abs(x))) - 0.5)

# Do simulation
for(m in 1:length(M)){
  for(s in 1:length(sigma)){
    t <- seq(0, 1, length.out = M[m] + 1) # Design points
    
    simulation_ <- rerun(N, quadratic_trajectory(M[m], sigma[s], hetero_noise, hetero_func))
    mean_ <- true_mean(t)
    covariance_ <- true_covariance(t, t)
    
    quadratic.trajectories <- list(
      simulation = simulation_,
      mean = mean_,
      covariance = covariance_
    )
    # Naming convention (fraction.brownian.trajectories-M-H-sigma)
    saveRDS(quadratic.trajectories, 
            file = paste0('./data/quadratic.trajectories-', M[m], '-', sigma[s], '-hetero3.rds'))
  }
}
