% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_b.R
\name{estimate_b}
\alias{estimate_b}
\title{Perform an estimation of the bandwidth for the smoothing of curves.}
\usage{
estimate_b(data, H0 = 0.5, L0 = 1, sigma = 0, K = "epanechnikov")
}
\arguments{
\item{data}{A list, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} The sampling points
\item \strong{$x} The observed points.
}}

\item{H0}{Numeric, an estimation of \eqn{H_0}.}

\item{L0}{Numeric, an estimation of \eqn{L_0}.}

\item{sigma}{Numeric, an estimation of \eqn{\sigma}.}

\item{K}{Character string, the kernel used for the estimation:
\itemize{
\item epanechnikov (default)
\item beta
\item uniform
}}
}
\value{
Numeric, an estimation of the bandwidth.
}
\description{
This function performs an estimation of the bandwidth for a univariate kernel
regression estimator defined over continuous data using the method of
\cite{add ref}. An estimation of \eqn{H_0}, \eqn{L_0} and \eqn{\sigma} have
to be provided to estimate the bandwidth.
}
\examples{
X <- generate_fractional_brownian(N = 1000, M = 300, H = 0.5, sigma = 0.05)
estimate_b(X, H0 = 0.5, L0 = 1, sigma = 0.05)
}
\seealso{
Other estimate bandwidth: 
\code{\link{estimate_b_cv}()},
\code{\link{estimate_b_list}()}
}
\concept{estimate bandwidth}
