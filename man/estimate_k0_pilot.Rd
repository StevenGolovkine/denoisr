% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_k0.R
\name{estimate_k0_pilot}
\alias{estimate_k0_pilot}
\title{Perform the estimation of the pilot \eqn{k_0}}
\usage{
estimate_k0_pilot(M)
}
\arguments{
\item{M}{Numeric, mean number of sampling points per curve}
}
\value{
Numeric, the pilot \eqn{k_0}
}
\description{
This function performs the estimation of the pilot \eqn{k_0} as used in
\cite{add ref}.
}
\examples{
estimate_k0_pilot(200)
}
\seealso{
Other estimate \eqn{k_0}: 
\code{\link{estimate_k0_oracle}()},
\code{\link{estimate_k0}()}
}
\concept{estimate \eqn{k_0}}
