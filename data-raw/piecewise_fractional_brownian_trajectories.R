######################################################################################
#                  Generate piecewise fractional Brownian motion                     #
######################################################################################

# Define functions

#' Generate piecewise fractional Brownian motion with a random noise.
#' 
#' @param M Number of points in the trajectory
#' @param H Vector of Hurst coefficients
#' @param sigma Standard deviation of the noise to add to the trajectory
#' @return A tibble containing the trajectory and the sampling points.
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
    x[i] <- list(as.vector(fbm(hurst = H[i], n = M_nn[i] - 1)))
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
  
  return(as_tibble(simu, .name_repair = 'unique'))
}


# Define some parameters
N <- 500000 # Number of curves
M <- c(1000)  # Number of points per curves (do it with 50, 200, 1000)
H <- c(0.2, 0.5, 0.8) # Hurst coefficient (do it with 0.4, 0.5, 0.6)
sigma <- c(0.01, 0.05, 0.1) # Standard deviation of the noise (do it with 0, 0.01, 0.05, 0.1)

# Do simulation
for(m in 1:length(M)){

      simulation_ <- rerun(N, piecewise_fractional_brownian_trajectory(M[m], H, sigma))
      
      fractional.brownian.trajectories <- list(
        simulation = simulation_
      )
      
      # Naming convention (fraction.brownian.trajectories-M-H-sigma)
      saveRDS(fractional.brownian.trajectories, 
              file = paste0('./data/piecewise.fractional.brownian.trajectories-', 
                            M[m], '-', paste(H, collapse = '.'), '-', paste(sigma, collapse = '.'), '.rds'))
}

