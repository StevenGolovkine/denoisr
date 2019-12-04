#' ######################################################################################
#' #       Simulation from Uniform Convergence Rates for Nonparametric Regression       #
#' #         and Principal Components Analysis in Functional / Longitudinal Data        #
#' #                                      Li and Hsing                                  #
#' ######################################################################################
#' 
#' # Load packages
#' library(tidyverse)
#' 
#' # Define functions
#' 
#' #' Simulation 1 from Li and Hsing.
#' #' 
#' #' @param M Number of design points
#' #' @param omega Vector of size 3 containing the variance of the factor
#' #' @param sigma Variance of the noise
#' #' @return A tibble containing the trajectories and the sampling points.
#' simulation_1 <- function(M, omega, sigma=0){
#'   # Grid of design points
#'   grid <- runif(M, 0, 1)
#'   grid <- grid[order(grid)]
#'   
#'   # Generate scores
#'   xi_1 <- rnorm(1, 0, sqrt(omega[1]))
#'   xi_2 <- rnorm(1, 0, sqrt(omega[2]))
#'   xi_3 <- rnorm(1, 0, sqrt(omega[3]))
#'   
#'   # Generate some noise
#'   epsilon <- rnorm(M, 0, sigma)
#'   
#'   # Generate the mean function and the eigenfunctions
#'   mu <- 5 * (grid - 0.6) ** 2
#'   phi_1 <- rep(1, M)
#'   phi_2 <- sqrt(2) * sin(2*pi*grid)
#'   phi_3 <- sqrt(2) * cos(2*pi*grid)
#'   
#'   # Generate the trajectory
#'   Y <- mu + xi_1 * phi_1 + xi_2 * phi_2 + xi_3 * phi_3 + epsilon
#'   
#'   return(tibble(t = grid, x = Y))
#' }
#' 
#' #' Simulation 2 from Li and Hsing.
#' #' 
#' #' @param M Number of design points
#' #' @param k Number of eigenelements considered
#' #' @param sigma Variance of the noise
#' #' @return A tibble containing the trajectories and the sampling points.
#' simulation_2 <- function(M, k, sigma){
#'   # Grid of design points
#'   grid <- runif(M, 0, 1)
#'   grid <- grid[order(grid)]
#'   
#'   # Generate some noise
#'   epsilon <- rnorm(M, 0, sigma)
#'   
#'   # Compute the trajectory
#'   Y <- rep(0, M)
#'   for(i in 1:k){
#'     w <- 4 / ((2*i - 1) ** 2 * pi ** 2)
#'     xi <- rnorm(1, 0, sqrt(w))
#'     phi <- sqrt(2) * sin((i - 0.5)*pi*grid)
#'     Y <- Y + w * phi
#'   }
#'   Y <- Y + epsilon
#'   
#'   return(tibble(t = grid, x = Y))
#' }
#' 
#' #' Compute the true mean of the Li trajectories for simulation 1.
#' #' 
#' #' @param t A vector of sampling points
#' #' @return A tibble containing the true mean trajectory and the sampling points. 
#' true_mean_1 <- function(t){
#'   x <- 5 * (t - 0.6) ** 2
#'   return(tibble(t = t, x = x))
#' }
#' 
#' #' Compute the true mean of the Li trajectories for simulation 2.
#' #' 
#' #' @param t A vector of sampling points
#' #' @return A tibble containing the true mean trajectory and the sampling points. 
#' true_mean_2 <- function(t){
#'   x <- rep(0, length(t))
#'   return(tibble(t = t, x = x))
#' }
#' 
#' #' Compute the true covariance of the Li trajectories for simulation 1.
#' #' 
#' #' @param s_ A vector of sampling points
#' #' @param t_ A vector of sampling points
#' #' @param omega Vector of size 3 containing the variance of the factor
#' #' @return A tibble containing the true covariance for each pair (s, t)
#' true_covariance_1 <- function(s_, t_, omega){
#'   phi <- matrix(rep(0, times = length(s_)*length(t_)), ncol = length(s_))
#'   for(is in seq_along(s_)){
#'     for(it in seq_along(t_)){
#'       phi[is, it] <- omega[1] + 2 * omega[2] * sin(2*pi*s_[is]) * sin(2*pi*t_[it]) +
#'                      2 * omega[3] * cos(2*pi*s_[is]) * cos(2*pi*t_[it])
#'     }
#'   }
#'   return(tibble(
#'     s = rep(s_, times = length(t_)),
#'     t = rep(t_, each = length(s_)),
#'     phi = as.vector(phi)))
#' }
#' 
#' 
#' #' Compute the true covariance of the Li trajectories for simulation 1.
#' #' 
#' #' @param s_ A vector of sampling points
#' #' @param t_ A vector of sampling points
#' #' @return A tibble containing the true covariance for each pair (s, t)
#' true_covariance_2 <- function(s_, t_){
#'   res <- tibble(
#'     s = rep(s_, times = length(t_)),
#'     t = rep(t_, each = length(s_))) %>% 
#'       mutate(phi = pmin(s, t))
#'   return(res)
#' }
#' 
#' # Define some parameters
#' N <- 100000 # Number of curves
#' M <- 1000 # Number of points per curves
#' omega <- c(0.6, 0.3, 0.1) # Variance of the eigenvalues
#' sigma <- 0.1 # Standard deviation of the noise
#' k <- 3 # Number of eigenelements to consider
#' t <- seq(0, 1, length.out = M+1) # Design points
#' 
#' # Do simulation
#' # Simulation 1
#' simulation1 <- rerun(N, simulation_1(M, omega, sigma))
#' mean_1 <- true_mean_1(t)
#' covariance_1 <- true_covariance_1(t, t, omega)
#' 
#' quadratic.trajectories.1000.0.1 <- list(
#'   simulation = simulation1,
#'   mean = mean_1,
#'   covariance = covariance_1
#' )
#' 
#' # Simulation 2
#' simulation2 <- rerun(N, simulation_2(M, k, sigma))
#' mean_2 <- true_mean_2(t)
#' covariance_2 <- true_covariance_2(t, t)
#' 
#' li.trajectories.2 <- list(
#'   simulation = simulation2,
#'   mean = mean_2,
#'   covariance = covariance_2
#' )
#' 
#' # Save data
#' usethis::use_data(quadratic.trajectories.1000.0.1, overwrite = TRUE)
#' usethis::use_data(li.trajectories.2, overwrite = TRUE)
