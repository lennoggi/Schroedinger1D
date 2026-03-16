#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <cmath>


// ***** DON'T TOUCH *****
#define GAUSSIAN 0
#define BOX      1
// ***********************

// Wave function type
#define WF GAUSSIAN


// Reduced Planck constant
constexpr inline double HBAR = 1.;

// Mass of the particle
constexpr inline double M = 1.;


// Number of spatial grid points
constexpr inline size_t NX = 1000;

// Spatial bounds are [0, L - L/N]
constexpr inline double L = 1.;


// Number of time steps
constexpr inline size_t NT = 100;

// Time step
constexpr inline double DT = 1.e-06;


// Initial wavefunction center
constexpr inline double X0 = 0.1*L;

// Initial wavefunction half-width
constexpr inline double SIGMA = 0.03*L;

// ***** DON'T TOUCH *****
// Half the maximum momentum of the wave function
constexpr inline double PMAX_HALF = M_PI*HBAR*static_cast<double>(NX)/L;
static_assert(PMAX_HALF > 0.);
// ***********************

// Initial wavefunction momentum
constexpr inline double P0 = 0.05*PMAX_HALF;


// ***** DON'T TOUCH *****
constexpr inline int FREE_PROPAGATION = 0;
constexpr inline int BARRIER_WELL     = 1;
constexpr inline int STEP             = 2;
constexpr inline int HARMONIC         = 3;

// Maximum allowed value for the potential
constexpr inline double VMAX = 4.*M_PI*HBAR/DT;
static_assert(VMAX > 0.);
// ***********************

// Potential choice
constexpr inline int POT = FREE_PROPAGATION;

// Central value for the potential
constexpr inline double VC = 0.5*L;

// Height of the potential barrier/well or step
constexpr inline double VH = 0.1*VMAX;

// Width of the potential barrier/well
constexpr inline double VW = 0.1*L;

// Strength of the harmonic potential
constexpr inline double K = 1.e+05;


// Output frequency
constexpr inline size_t OUT_EVERY = 5;

// Output filename
#define FILENAME "Data.h5"

// Verbosity
constexpr inline bool VERBOSE = true;


#endif  // PARAMETERS_HH
