#include "include/Declare_functions.hh"
#include "Parameters.hh"

double step(const double &x) {
    return (x >= VC) ? VH : 0.;
}

double barrier_well(const double &x) {
    constexpr double x1 = VC - VW/2.;
    constexpr double x2 = VC + VW/2.;
    return (x >= x1 and x <= x2) ? VH : 0.;
}

double harmonic(const double &x) {
    constexpr double strength = K/2.;
    const     double xs       = x - VC;
    return strength*xs*xs;
}
