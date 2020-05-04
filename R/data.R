################################################################################
#                              Document the data                               #           
################################################################################

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