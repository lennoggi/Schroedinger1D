#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <complex>

std::complex<double> gaussian_wf (const double &x);
std::complex<double> gaussian_Fwf(const double &p);
std::complex<double>      box_wf (const double &x);
std::complex<double>      box_Fwf(const double &p);

#endif  // DECLARE_FUNCTIONS_HH
