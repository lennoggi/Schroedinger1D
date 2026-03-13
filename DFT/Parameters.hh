#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <cmath>

// ***** DON'T TOUCH *****
#define GAUSSIAN 0
#define BOX      1
// ***********************

// Wave function type
#define WF BOX

// Reduced Planck constant
constexpr inline double HBAR = 1.;

// Number of spatial grid points
constexpr inline int N = 1000;

// Spatial bounds are [0, L - L/N]
constexpr inline double L = 1.;

// Half the maximum momentum of the wave function
constexpr inline double PMAX_HALF = M_PI*HBAR*static_cast<double>(N)/L;

// Initial wavefunction center
constexpr inline double X0 = 0.1*L;

// Initial wavefunction half-width
constexpr inline double SIGMA = 0.03*L;

// Initial wavefunction momentum
constexpr inline double P0 = 0.05*PMAX_HALF;

// Output filename
#define FILENAME "Data.h5"

#endif  // PARAMETERS_HH
