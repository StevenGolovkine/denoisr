######################################################################################
#               Simulation from Common Functional Principal Components               #
#                               Benko, Härdle and Kneip                              #
######################################################################################

# Load packages
library(tidyverse)

# Define functions

#' Simulation named 'setup a' from Benko, Härdle and Kneip.
#' 
#' @param M Number of design points
#' @param lambda Vector of size 4 containing the variance of the factor
#' @param delta Shift parameter between the function
#' @return A tibble containing the trajectories and the sampling points.
setup_a <- function(M, lambda, delta){
  # Grid of design points
  t <- seq(1, M, by = 1) / M 
  
  # Generate factor loadings
  beta_11 <- rnorm(1, mean = 0, sd = sqrt(lambda[1]))
  beta_21 <- rnorm(1, mean = 0, sd = sqrt(lambda[2]))
  beta_12 <- rnorm(1, mean = 0, sd = sqrt(lambda[3]))
  beta_22 <- rnorm(1, mean = 0, sd = sqrt(lambda[4]))
  
  # Generate the trajectories
  X_1 <- sqrt(2) * (beta_11 * sin(2*pi*t) + beta_12 * cos(2*pi*t))
  X_2 <- sqrt(2) * (beta_21 * sin(2*pi*(t + delta)) + beta_22 * cos(2*pi*(t + delta)))
  
  return(tibble(t = t, x = X_1, y = X_2))
}

#' Simulation named 'setup b' from Benko, Härdle and Kneip.
#' 
#' @param M Number of design points
#' @param lambda Vector of size 4 containing the variance of the factor
#' @param delta Shift parameter between the function
#' @return A tibble containing the trajectories and the sampling points.
setup_b <- function(M, lambda, delta){
  # Grid of design points
  t <- seq(1, M, by = 1) / M 
  
  # Generate factor loadings
  beta_11 <- rnorm(1, mean = 0, sd = sqrt(lambda[1]))
  beta_21 <- rnorm(1, mean = 0, sd = sqrt(lambda[2]))
  beta_12 <- rnorm(1, mean = 0, sd = sqrt(lambda[3]))
  beta_22 <- rnorm(1, mean = 0, sd = sqrt(lambda[4]))
  
  # Generate the trajectories
  X_1 <- sqrt(2) * (beta_11 * sin(2*pi*t) + beta_12 * cos(2*pi*t))
  X_2 <- sqrt(2) * (beta_21 * sin(2*pi*(t + delta)) + beta_22 * sin(4*pi*(t + delta)))
  
  return(tibble(t = t, x = X_1, y = X_2))
}

#' Simulation named 'setup noisy a' from Benko, Härdle and Kneip.
#' 
#' @param M Number of design points
#' @param lambda Vector of size 4 containing the variance of the factor
#' @param delta Shift parameter between the function
#' @param sigma Standard deviation of the noise
#' @return A tibble containing the trajectories and the sampling points.
setup_a_noise <- function(M, lambda, delta, sigma){
  # Grid of design points
  t <- seq(1, M, by = 1) / M 
  
  # Generate factor loadings
  beta_11 <- rnorm(1, mean = 0, sd = sqrt(lambda[1]))
  beta_21 <- rnorm(1, mean = 0, sd = sqrt(lambda[2]))
  beta_12 <- rnorm(1, mean = 0, sd = sqrt(lambda[3]))
  beta_22 <- rnorm(1, mean = 0, sd = sqrt(lambda[4]))
  
  # Generate some noise
  epsilon <- rnorm(length(t), 0, sigma)
  
  # Generate the trajectories
  X_1 <- sqrt(2) * (beta_11 * sin(2*pi*t) + beta_12 * cos(2*pi*t))
  X_2 <- sqrt(2) * (beta_21 * sin(2*pi*(t + delta)) + beta_22 * cos(2*pi*(t + delta)))
  
  return(tibble(t = t, x = X_1 + epsilon, y = X_2 + epsilon))
}

#' Compute the true mean of the Kneip trajectories for setup_a and setup_a_noise.
#' 
#' @param t A vector of sampling points
#' @return A tibble containing the true mean trajectory and the sampling points. 
true_mean_a <- function(t){
  x <- rep(0, length(t))
  y <- rep(0, length(t))
  return(tibble(t = t, x = x, y = y))
}

#' Compute the true mean of the Kneip trajectories for setup_b.
#' 
#' @param t A vector of sampling points
#' @return A tibble containing the true mean trajectory and the sampling points. 
true_mean_b <- function(t){
  x <- rep(0, length(t))
  y <- rep(0, length(t))
  return(tibble(t = t, x = x, y = y))
}

#' Compute the true covariance of the Kneip trajectories for setup_a and
#' setup_a_noise.
#' 
#' @param s_ A vector of sampling points
#' @param t_ A vector of sampling points
#' @return A tibble containing the true covariance for each pair (s, t)
true_covariance_a <- function(s_, t_, lambda, delta){
  phi_x <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
  phi_y <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
  
  for(is in seq_along(s_)){
    for(it in seq_along(t_)){
      phi_x[is, it] <- 2 * lambda[1] * sin(2*pi*s_[is]) * sin(2*pi*t_[it]) + 
                       2 * lambda[2] * cos(2*pi*s_[is]) * cos(2*pi*t_[it])
      phi_y[is, it] <- 2 * lambda[3] * sin(2*pi*(s_[is] + delta)) * sin(2*pi*(t_[it] + delta)) + 
                       2 * lambda[4] * cos(2*pi*(s_[is] + delta)) * cos(2*pi*(t_[it] + delta))
    }
  }
  return(tibble(
    s = rep(s_, times = length(t_)),
    t = rep(t_, each = length(s_)),
    phi_x = as.vector(phi_x),
    phi_y = as.vector(phi_y)))
}

#' Compute the true covariance of the Kneip trajectories for setup_b.
#' 
#' @param s_ A vector of sampling points
#' @param t_ A vector of sampling points
#' @return A tibble containing the true covariance for each pair (s, t)
true_covariance_b <- function(s_, t_, lambda, delta){
  phi_x <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
  phi_y <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
  
  for(is in seq_along(s_)){
    for(it in seq_along(t_)){
      phi_x[is, it] <- 2 * lambda[1] * sin(2*pi*s_[is]) * sin(2*pi*t_[it]) + 
                       2 * lambda[2] * cos(2*pi*s_[is]) * cos(2*pi*t_[it])
      phi_y[is, it] <- 2 * lambda[3] * sin(2*pi*(s_[is] + delta)) * sin(2*pi*(t_[it] + delta)) + 
                       2 * lambda[4] * sin(4*pi*(s_[is] + delta)) * sin(4*pi*(t_[it] + delta))
    }
  }
  return(tibble(
    s = rep(s_, times = length(t_)),
    t = rep(t_, each = length(s_)),
    phi_x = as.vector(phi_x),
    phi_y = as.vector(phi_y)))
}

# Define some parameters
N <- 100 # Number of curves
M <- 50 # Number of points per curves
lambda <- c(10, 5, 8, 4) # Variance of the eigenvalues
delta <- 0.05 # Shift
sigma <- sqrt(0.25) # Standard deviation of the noise
t <- seq(0, 1, length.out = M+1) # Design points

# Do simulation
# Setup (a)
simulation_a <- rerun(N, setup_a(M, lambda, delta))
mean_a <- true_mean_a(t)
covariance_a <- true_covariance_a(t, t, lambda, delta)

kneip.trajectories.a <- list(
  simulation = simulation_a,
  mean = mean_a,
  covariance = covariance_a
)

# Setup (a_noise)
simulation_a_noise <- rerun(N, setup_a_noise(M, lambda, delta, sigma))
mean_a_noise <- true_mean_a(t)
covariance_a_noise <- true_covariance_a(t, t, lambda, delta)

kneip.trajectories.a.noise <- list(
  simulation = simulation_a_noise,
  mean = mean_a_noise,
  covariance = covariance_a_noise
)

# Setup (b)
simulation_b <- rerun(N, setup_b(M, lambda, delta))
mean_b <- true_mean_b(t)
covariance_b <- true_covariance_b(t, t, lambda, delta)

kneip.trajectories.b <- list(
  simulation = simulation_b,
  mean = mean_b,
  covariance = covariance_b
)

# Save data
usethis::use_data(kneip.trajectories.a, overwrite = TRUE)
usethis::use_data(kneip.trajectories.a.noise, overwrite = TRUE)
usethis::use_data(kneip.trajectories.b, overwrite = TRUE)

