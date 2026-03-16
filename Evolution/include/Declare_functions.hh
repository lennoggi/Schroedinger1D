#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <complex>

// Wave function and their Fourier transforms
std::complex<double> gaussian_wf (const double &x);
std::complex<double> gaussian_Fwf(const double &p);
std::complex<double>      box_wf (const double &x);
std::complex<double>      box_Fwf(const double &p);

// Potentials
double barrier_well(const double &x);
double         step(const double &x);
double     harmonic(const double &x);

#endif  // DECLARE_FUNCTIONS_HH
