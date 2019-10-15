######################################################################################
#         Simulation from From Sparse to Dense Functional Data and Beyond            #
#                                   Zhang and Wang                                   #
######################################################################################

# Load packages
library(tidyverse)

# Define functions

#' Simulation from Zhang and Wang.
#' 
#' @param M Number of design points
#' @param sigma Variance of the noise
#' @return A tibble containing the trajectories and the sampling points.
zhang_trajectory <- function(M, sigma){
  # Grid of design points
  grid <- runif(M, 0, 1)
  grid <- c(0, grid[order(grid)])
  
  # Generate scores
  A_1 <- rnorm(1, 0, sqrt(1/4))
  A_2 <- rnorm(1, 0, sqrt(1/9))
  A_3 <- rnorm(1, 0, sqrt(1/16))
  A_4 <- rnorm(1, 0, sqrt(1/25))
  
  # Generate some noise
  epsilon <- rnorm(M+1, 0, sigma)
  
  # Compute the mean function and the eigenfunctions
  mu <- 1.5 * sin(3*pi*(grid + 0.5)) + 2 * grid ** 3
  phi_1 <- sqrt(2) * cos(2*pi*grid)
  phi_2 <- sqrt(2) * sin(2*pi*grid)
  phi_3 <- sqrt(2) * cos(4*pi*grid)
  phi_4 <- sqrt(2) * sin(4*pi*grid)
  
  # Generate the trajectory
  Y <- mu + A_1 * phi_1 + A_2 * phi_2 + A_3 * phi_3 + A_4 * phi_4 + epsilon
  
  return(tibble(t = grid, x = Y))
}

#' Simulate the number of grid points.
#' 
#' @param N Number of curves to simulate
#' @param setting The setting to follow
#' @return The number of grid points for a curve
simulation_M <- function(N, setting = 1){
  if(setting == 1){
    M <- sample(c(2, 3, 4, 5), 1)
  } else if(setting == 2){
    M <- sample(c(N/4, sample(c(2, 3, 4, 5), 1)), 1, prob = c(0.5, 0.5))
  } else if(setting == 3){
    M <- sample(c(N/4, sample(c(2, 3, 4, 5), 1)), 1, prob = c(1/sqrt(N), 1 - 1/sqrt(N)))
  } else if(setting == 4){
    M <- sample(seq(round(N/8,0), round(3*N/8,0), 1), 1)
  } else{
    cat('Setup does not exist!')
    M <- 1
  }
  return(M)
}

#' Compute the true mean of the simulation from Zhang and Wang.
#' 
#' @param t A vector of sampling points
#' @return A tibble containing the true mean trajectory and the sampling points. 
true_mean <- function(t){
  res <- mu <- 1.5 * sin(3*pi*(t + 0.5)) + 2 * t ** 3
  return(tibble(t = t, x = res))
}

#' Compute the true covariance of the simulation from Zhang and Wang.
#' 
#' @param s_ A vector of sampling points
#' @param t_ A vector of sampling points
#' @return A tibble containing the true covariance for each pair (s, t)
true_covariance <- function(s_, t_){
  phi_1 <- sqrt(2) * cos(2*pi*s_)
  phi_2 <- sqrt(2) * sin(2*pi*s_)
  phi_3 <- sqrt(2) * cos(4*pi*s_)
  phi_4 <- sqrt(2) * sin(4*pi*s_)
  
  phi_1t <- t(sqrt(2) * cos(2*pi*t_))
  phi_2t <- t(sqrt(2) * sin(2*pi*t_))
  phi_3t <- t(sqrt(2) * cos(4*pi*t_))
  phi_4t <- t(sqrt(2) * sin(4*pi*t_))
  
  phi <- (phi_1 %*% phi_1t) / 4 + (phi_2 %*% phi_2t) / 9 + 
         (phi_3 %*% phi_3t) / 16 + (phi_4 %*% phi_4t) / 25 
  return(tibble(
    s = rep(s_, times = length(t_)),
    t = rep(t_, each = length(s_)),
    phi = as.vector(phi)))
}


# Define some parameters
N <- 10000 # Number of curves
sigma <- sqrt(0.1) # Standard deviation of the noise
t <- seq(0, 1, length.out = 100) # Design points

# Do simulation
# Setting 1
simulation <- rerun(N, {simulation_M(N, setting = 1) %>% zhang_trajectory(sigma)})
mean_ <- true_mean(t)
covariance_ <- true_covariance(t, t)

zhang.trajectories.1 <- list(
  simulation = simulation,
  mean = mean_,
  covariance = covariance_
)

# Setting 2
simulation <- rerun(N, {simulation_M(N, setting = 2) %>% zhang_trajectory(sigma)})
mean_ <- true_mean(t)
covariance_ <- true_covariance(t, t)

zhang.trajectories.2 <- list(
  simulation = simulation,
  mean = mean_,
  covariance = covariance_
)

# Setting 3
simulation <- rerun(N, {simulation_M(N, setting = 3) %>% zhang_trajectory(sigma)})
mean_ <- true_mean(t)
covariance_ <- true_covariance(t, t)

zhang.trajectories.3 <- list(
  simulation = simulation,
  mean = mean_,
  covariance = covariance_
)

# Setting 4
simulation <- rerun(N, {simulation_M(N, setting = 4) %>% zhang_trajectory(sigma)})
mean_ <- true_mean(t)
covariance_ <- true_covariance(t, t)

zhang.trajectories.4 <- list(
  simulation = simulation,
  mean = mean_,
  covariance = covariance_
)


# Save data
usethis::use_data(zhang.trajectories.1, overwrite = TRUE)
usethis::use_data(zhang.trajectories.2, overwrite = TRUE)
usethis::use_data(zhang.trajectories.3, overwrite = TRUE)
usethis::use_data(zhang.trajectories.4, overwrite = TRUE)