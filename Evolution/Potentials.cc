#include "include/Declare_functions.hh"
#include "Parameters.hh"

double step(const double &x) {
    return (x >= VC) ? VH : 0.0;
}

double barrier_well(const double &x) {
    constexpr double vw_half = VW/2.0;
    constexpr double x1 = VC - vw_half;
    constexpr double x2 = VC + vw_half;
    return (x >= x1 and x <= x2) ? VH : 0.0;
}

double harmonic(const double &x) {
    constexpr double strength = K/2.0;
    const     double xs       = x - VC;
    return strength*xs*xs;
}
