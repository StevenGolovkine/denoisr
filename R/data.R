################################################################################
#                              Document the data                               #           
################################################################################

#' Fractional brownian trajectories
#' 
#' A dataset containing the simulation of fractional Brownian motions. The 
#' increments of a fractional Brownian motion are not independent. A fractional 
#' Brownian motion is characterized by a parameter \eqn{H}, named Hurst 
#' coefficient. They are defined on \eqn{[0, 1]}.
#' 
#' @format A list with 1000 elements. Each element corresponds to one simulation of 
#' a fractional Brownian trajectory. Each element is a list with 3 elements which
#' are:
#' \itemize{
#'  \item \strong{t} The sampling points
#'  \item \strong{x} The trajectory with random noise
#'  \item \strong{x_true} The true underlying trajectory
#' }
"fractional_brownian"

#' Piecewise fractional brownian trajectories
#' 
#' A dataset containing the simulation of piecewise fractional Brownian motions.
#' A piecewise fractional Brownian motion is defined by a non constant Hurst 
#' parameter along the sampling points. We observe the process at regularly 
#' spaced time \eqn{t_i = \frac{i}{M_n}}, where \eqn{i = 0, \dots, M_n}.
#' We define a segmentation \eqn{\tau = (\tau_k)_{k=0, \dots, K+1}}, with 
#' \eqn{0 = \tau_0 < \tau_1 < \dots < \tau_{K} < \tau_{K+1} = 1}. So, on the 
#' interval \eqn{[\tau_k, \tau_{k+1}]}, for \eqn{k = 0, \dots, K}, the process 
#' is a fractional Brownian motion with Hurst parameter \eqn{H_k}. 
#' 
#' @format A list with 1000 elements. Each element corresponds to one simulation of 
#' a piecewise fractional Brownian trajectory. Each element is a list with 
#' 3 elements which are:
#' \itemize{
#'  \item \strong{t} The sampling points
#'  \item \strong{x} The trajectory with random noise
#'  \item \strong{x_true} The true underlying trajectory
#' }
"piecewise_fractional_brownian"

#' Canadian average annual temperature
#' 
#' Daily temperature at 35 different locations in Canada averaged over 1960 to
#' 1994.
#' 
#' @format A list with 35 elements. Each element corresponds to a particular 
#' Canadian station. Each station is a list with 2 elements which are:
#' \itemize{
#'  \item \strong{t} The day of the year the temperature is taken (normalized on \eqn{[0, 1]})
#'  \item \strong{x} The average temperature for each day of the year
#' }
#'
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.
#' @references Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer, New York
#' 
#' @seealso \code{\link[fda]{CanadianWeather}}
"canadian_temperature_daily"