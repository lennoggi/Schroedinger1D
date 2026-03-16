#include <cmath>
#include "include/Declare_functions.hh"
#include "Parameters.hh"

double barrier_well(const double &x) {
    constexpr double half_width = VC - VW/2.;
    return (fabs(x) <= half_width) ? VH : 0.;
}

double step(const double &x) {
    return (x >= VC) ? VH : 0.;
}

double harmonic(const double &x) {
    constexpr double strength = K/2.;
    const     double xs       = x - VC;
    return strength*xs*xs;
}
