#include <cstdlib>
#include <cmath>
#include <complex>

#include "include/Declare_functions.hh"
#include "Parameters.hh"

using namespace std;
using namespace std::complex_literals;


// Gaussian
complex<double> gaussian_wf(const double &x) {
    constexpr double norm         = pow(2.0/(M_PI*SIGMA*SIGMA), 0.25);
    constexpr double p0_over_hbar = P0/HBAR;
    const     double           xs = x - X0;
    const     double           xx = xs/SIGMA;
    return norm*exp(-xx*xx - 1.0i*p0_over_hbar*xs);
}

complex<double> gaussian_Fwf(const double &p) {
    constexpr double norm         = pow(SIGMA*SIGMA/(2.0*M_PI*HBAR*HBAR), 0.25);
    constexpr double x0_over_hbar = X0/HBAR;
    constexpr double sigmabar     = SIGMA/(2.0*HBAR);
    const     double pp           = sigmabar*(p + P0);
    return norm*exp(-pp*pp - 1.0i*x0_over_hbar*p);
}


// Box
complex<double> box_wf(const double &x) {
    constexpr double norm         = 1.0/sqrt(2.0*SIGMA);
    constexpr double p0_over_hbar = P0/HBAR;
    const     double xs           = x - X0;
    return (fabs(xs) <= SIGMA) ? norm*exp(-1.0i*p0_over_hbar*xs) : 0.0;
}

complex<double> box_Fwf(const double &p) {
    constexpr double sigma_over_hbar = SIGMA/HBAR;
    constexpr double norm            = 1.0/sqrt(M_PI*sigma_over_hbar);
    constexpr double x0_over_hbar    = X0/HBAR;
    const     double ps              = p - P0;
    return norm*exp(-1.0i*x0_over_hbar*p)*sin(sigma_over_hbar*ps)/ps;
}
