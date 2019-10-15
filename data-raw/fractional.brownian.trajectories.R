######################################################################################
#                       Generate fractional Brownian motion                          #
######################################################################################

# Load packages
library(tidyverse)

# Define functions

#' Generate fractional Brownian motion with a random noise.
#' 
#' @param M Number of points in the trajectory
#' @param sigma Standard deviation of the noise to add to the trajectory, can be a list.
#' @param H Hurst coefficient
#' @param hetero_func A function to generate the noise
#' @return A tibble containing the trajectory and the sampling points.
fractional_brownian_trajectory <- function(M, H, sigma, hetero_func = NULL){
  require(somebm)
  if(is.null(hetero_func)){
    sd.function <- function(t) (abs(t)+1)**(5/6)
  } else{
    sd.function <- hetero_func
  }
  
  M_n <- rpois(1, M)
  t <- seq(0, 1, length.out = M_n + 1)
  x <- as.vector(fbm(hurst = H, n = M_n))
  
  # Start to fill the data
  simu <- matrix(rep(0, (M_n + 1) * (length(sigma) + 3)), nrow = M_n + 1)
  simu[, 1] <- t
  simu[, 2] <- x
  
  e <- rnorm(M_n + 1, mean = 0, sd = 1)
  
  # Add columns with homoscedastic noise.
  j <- 3
  for(i in sigma){
    simu[, j] <- x + i * e
    j = j + 1
  }
  
  # Add a column with heteroscedastique noise.
  e <- sd.function((x - median(x)) / sd(x)) * rnorm(M_n + 1, mean = 0, sd = 0.05)
  simu[, ncol(simu)] <- x + e
  
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
N <- 10000 # Number of curves
M <- c(50, 200, 1000)  # Number of points per curves (do it with 50, 200, 1000)
H <- c(0.9) # Hurst coefficient (do it with 0.4, 0.5, 0.6)
sigma <- c(0.01, 0.05, 0.1, 0.25) # Standard deviation of the noise (do it with 0.01, 0.05, 0.1, 0.25)
hetero_func <- function(x) return (2 * exp(abs(x)) / (1 + exp(abs(x))) - 0.5)

# Do simulation
for(m in 1:length(M)){
    t <- seq(0, 1, length.out = M[m] + 1) # Design points
    
    simulation_ <- rerun(N, fractional_brownian_trajectory(M[m], H, sigma, hetero_func))
    mean_ <- true_mean(t)
    covariance_ <- true_covariance(t, t, H)
    
    fractional.brownian.trajectories <- list(
      simulation = simulation_,
      mean = mean_,
      covariance = covariance_
    )
    # Naming convention (fraction.brownian.trajectories-M-H-sigma)
    saveRDS(fractional.brownian.trajectories, 
            file = paste0('./data/fractional.brownian.trajectories-', M[m], '-', H, '-', paste(sigma, collapse = '.'), '.rds'))
}


# Save data
# Naming convention (fraction.brownian.trajectories-M-H-sigma)
#usethis::use_data(fractional.brownian.trajectories, overwrite = TRUE)

