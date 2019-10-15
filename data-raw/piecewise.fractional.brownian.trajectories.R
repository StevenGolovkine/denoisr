######################################################################################
#                  Generate piecewise fractional Brownian motion                     #
######################################################################################

# Load packages
library(tidyverse)


# Define functions

#' Generate piecewise fractional Brownian motion with a random noise.
#' 
#' @param M Number of points in the trajectory
#' @param H Vector of Hurst coefficients
#' @param sigma Standard deviation of the noise to add to the trajectory
#' @param hetero_noise Should the moise be heteroscedastic? (default = FALSE)
#' @param hetero_func A function to generate the noise
#' @param random_cut Random cut for the change of regularity?
#' @return A tibble containing the trajectory and the sampling points.
piecewise_fractional_brownian_trajectory <- function(M, H, sigma, hetero_noise = FALSE, hetero_func = NULL, random_cut = FALSE){
  require(somebm)
  if(is.null(hetero_func)){
    sd.function <- function(t) (abs(t)+1)**(5/6)
  } else{
    sd.function <- hetero_func
  }
  
  if(random_cut){
    M_n <- rpois(1, M)
    if(M_n < 40) M_n <- 40
    
    M_nn <- vector(length = length(H))
    for(i in 1:(length(H) - 1)){
      M_nn[i] <- round(runif(1, min = (M_n / length(H)) - 5, max = (M_n / length(H)) + 5))
    }
    M_nn[length(H)] <- M_n - sum(M_nn)    
  } else {
    M_n <- M
    
    M_nn <- vector(length = length(H))
    for(i in 1:(length(H) - 1)){
      M_nn[i] <- round(M_n / length(H))
    }
    M_nn[length(H)] <- M_n - sum(M_nn)
  }

  t <- seq(0, 1, length.out = M_n)
  x <- c()
  for(i in seq_along(H)){
    x <- c(x, as.vector(fbm(hurst = H[i], n = M_nn[i] - 1)))
  }

  e <- rnorm(M_n, mean = 0, sd = sigma)
  if(hetero_noise) e <- sd.function((x - median(x)) / sd(x)) * e
  
  return(tibble(t = t, x = x + e))
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
H <- c(0.6, 0.5, 0.4) # Hurst coefficient (do it with 0.4, 0.5, 0.6)
sigma <- c(0, 0.01, 0.05, 0.1) # Standard deviation of the noise (do it with 0, 0.01, 0.05, 0.1)
hetero_noise <- FALSE # Add heteroscedastic noise?
hetero_func <- function(x) return (2 * exp(abs(x)) / (1 + exp(abs(x))) - 0.5)


# Do simulation
for(m in 1:length(M)){
  for(s in 1:length(sigma)){
      t <- seq(0, 1, length.out = M[m] + 1) # Design points
      
      simulation_ <- rerun(N, piecewise_fractional_brownian_trajectory(M[m], H, sigma[s], hetero_noise, hetero_func))
      mean_ <- true_mean(t)
      covariance_ <- true_covariance(t, t, mean(H))
      
      fractional.brownian.trajectories <- list(
        simulation = simulation_,
        mean = mean_,
        covariance = covariance_
      )
      # Naming convention (fraction.brownian.trajectories-M-H-sigma)
      saveRDS(fractional.brownian.trajectories, 
              file = paste0('./data/piecewise.fractional.brownian.trajectories-', M[m], '-', paste(H, collapse = '.'), '-', sigma[s], '-3.rds'))
  }
}

