# ####
# #### Load packages
# library(np)
# library(microbenchmark)
# library(somebm)
# 
# Rcpp::sourceCpp('./src/estimate_curve.cpp')
# 
# ####
# #### Generate data
# M <- 1000
# t <- seq(0, 1, length.out = M + 1)
# x <- as.vector(bm(n = M))
# y <- x + sqrt(0.05) * rnorm(M + 1)
# 
# ####
# #### Analysis with Uniform
# bw_uni <- npregbw(y ~ t,
#                   regtype = 'lc', 
#                   bwmethod = 'cv.ls',
#                   ckertype = 'uniform')$bw
# 
# fit_uni <- npreg(y ~ t, bws = bw_uni, ckertype = 'uniform')
# fit2_uni <- uniKernelSmoothingCurve(t, t, y, b = bw_uni)
# 
# ####
# #### Analysis with Epanechnikov
# bw_epa <- npregbw(y ~ t, 
#                   regtype = 'lc', 
#                   bwmethod = 'cv.ls', 
#                   ckertype = 'epanechnikov')$bw
# 
# # Careful sqrt(5) factor for the Epanechnikov kernel in banwitdth estimation
# fit_epa <- npreg(y ~ t, bws = bw_epa, ckertype = 'epanechnikov')
# fit2_epa <- epaKernelSmoothingCurve(t, t, y, b = sqrt(5)*bw_epa)
# 
# 
# ####
# #### Plots
# 
# # Uniform
# plot(t, y, col = 'black', cex = .5)
# lines(t, x, col = 'red')
# lines(t, fit_uni$mean, col = 'blue')
# lines(t, fit2_uni[,1], col = 'green')
# 
# # Epanechnikov
# plot(t, y, col = 'black', cex= .5)
# lines(t, x, col = 'red')
# lines(t, fit_epa$mean, col = 'blue')
# lines(t, fit2_epa[,1], col = 'green')
# 
# ####
# #### Test speed numerical integration
# f <- function(x, h){
#   return(0.75 * abs(1 - x**2) * abs(x)**h)
# }
# 
# h <- 0.7
# microbenchmark(
#   3 / ((h + 1) * (h + 3)),
#   integrate(f, -1, 1, h = h),
#   times = 1000
# )
# 
# ####
# #### Test speed npreg vs. mine
# microbenchmark(
#   npreg(y ~ t, bws = bw_uni, ckertype = 'uniform'),
#   uniKernelSmoothingCurve(t, t, y, b = bw_uni),
#   times = 100
# )

# # 
# # 
# # library(tictoc)
# # 
# # tic()
# # simulation2 <- simulation$simulation %>% 
# #   map(~ add_column(.x, 
# #                    ...6 = sort(runif(length(.x$...1), 0, 1)),
# #                    ...7 = sort(rbeta(length(.x$...1), 0.25, 0.25))))
# # toc()
# # 
# # simulation2 <- list(simulation2, simulation$mean, simulation$covariance, simulation$sigma)
# # 
# # saveRDS(simulation2, 
# #         file = paste0('./fractional.brownian.trajectories-1000-0.4-0.01.0.05.0.1-2.rds'))
