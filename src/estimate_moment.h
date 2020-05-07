// estimate_moment.h

#ifndef ESTIMATE_MOMENT_H
#define ESTIMATE_MOMENT_H

#include <RcppArmadillo.h>

arma::vec LOOmean(
    const List & curves, // Curves list ($x and $t)
    const arma::vec & U, // Estimation points
    const List & b, // Smoothing bandwidths
    const arma::uword & n // Curve to remove from the estimation
);

arma::mat covariance(
    const List & curves, // Curves list ($x and $t)
    const arma::vec & U, // Estimation points
    const List & b, // Smoothing bandwidth 
    const List & h // Smoothing bandwidth
);
  
#endif