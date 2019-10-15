// estimate_sigma.h

#ifndef ESTIMATE_SIGMA_H
#define ESTIMATE_SIGMA_H

#include <RcppArmadillo.h>

double estimateSigma(
    const List & curves // Curves list ($x and $t)
);

#endif