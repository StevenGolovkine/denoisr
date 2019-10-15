// estimtate_risk.h

#ifndef ESTIMATE_RISK_H
#define ESTIMATE_RISK_H

#include <RcppArmadillo.h>

List estimateRisk(
    const List & curves, // Curves list ($x and $t)
    const List & curves_estim, // Estimated curves list ($x and $t)
    const double & t0 // In which point estimate the risk.
);

List estimateRiskCurve(
    const List & curve, // Curve ($x and $t)
    const List & curve_estim // Estimated curve ($x and $t)
);

List estimateRiskCurves(
    const List & curves, // Curves list ($x and $t)
    const List & curves_estim // Estimated curves list ($x and $t)
);

#endif