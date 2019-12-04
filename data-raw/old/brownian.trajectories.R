#' ######################################################################################
#' #                         Generate standard Brownian motion                          #
#' ######################################################################################
#' 
#' # Load packages
#' library(tidyverse)
#' 
#' # Define functions
#' 
#' #' Generate standard Brownian motion with a random noise.
#' #' 
#' #' @param M Number of points in the trajectory
#' #' @param sigma Standard deviation of the noise to add to the trajectory
#' #' @param hetero_noise Should the moise be heteroscedastic? (default = FALSE)
#' #' @return A tibble containing the trajectory and the sampling points.
#' brownian_trajectory <- function(M, sigma, hetero_noise = FALSE, hetero_func = NULL){
#'   if(is.null(hetero_func)){
#'     sd.function <- function(t) (abs(t)+1)**(5/6)
#'   } else{
#'     sd.function <- hetero_func
#'   }
#'   
#'   M_n <- rpois(1, M)
#'   t <- c(0, sort(runif(M_n)))
#'   e <- rnorm(M_n + 1, mean = 0, sd = sigma)
#' 
#'   x <- cumsum(sapply(t - dplyr::lag(t, default = 0), FUN = function(s) rnorm(1, 0, s)))
#'   if(hetero_noise) e <- sd.function(x) * e
#'   return(tibble(t = t, x = x + e))
#' }
#' 
#' #' Compute the true mean of the generated Brownian motion.
#' #' 
#' #' @param t A vector of sampling points
#' #' @return A tibble containing the true mean trajectory and the sampling points. 
#' true_mean <- function(t){
#'   x <- rep(0, length(t))
#'   return(tibble(t = t, x = x))
#' }
#' 
#' #' Compute the true covariance of the generated Brownian motion.
#' #' 
#' #' @param s_ A vector of sampling points
#' #' @param t_ A vector of sampling points
#' #' @return A tibble containing the true covariance for each pair (s, t)
#' true_covariance <- function(s_, t_){
#'   res <- tibble(
#'     s = rep(s_, times = length(t_)),
#'     t = rep(t_, each = length(s_))) %>% 
#'     mutate(phi = pmin(s, t))
#'   return(res)
#' }
#' 
#' # Define some parameters
#' N <- 100000 # Number of curves
#' M <- 50 # Number of points per curves
#' sigma <- 0 # Standard deviation of the noise
#' t <- seq(0, 1, length.out = M) # Design points
#' hetero_func <- function(x) return (2 * exp(abs(x)) / (1 + exp(abs(x))) - 0.5)
#' 
#' # Do simulation
#' simulation_ <- rerun(N, brownian_trajectory(M, sigma))
#' mean_ <- true_mean(t)
#' covariance_ <- true_covariance(t, t)
#' 
#' brownian.trajectories <- list(
#'   simulation = simulation_,
#'   mean = mean_,
#'   covariance = covariance_
#' )
#' 
#' # Save data
#' usethis::use_data(brownian.trajectories, overwrite = TRUE)
