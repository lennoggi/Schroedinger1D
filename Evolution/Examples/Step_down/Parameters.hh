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
constexpr inline double HBAR = 1.0;

// Mass of the particle
constexpr inline double M = 1.0;


// Number of spatial grid points
constexpr inline size_t NX = 1000;

// Spatial bounds are [0, L - L/N]
constexpr inline double L = 1.0;


// Number of time steps
constexpr inline size_t NT = 200;

// Time step
constexpr inline double DT = 1.0e-06;


// Initial wavefunction center
constexpr inline double X0 = 0.45*L;

// Initial wavefunction half-width
constexpr inline double SIGMA = 0.03*L;

// ***** DON'T TOUCH *****
// Half the maximum momentum of the wave function
constexpr inline double PMAX_HALF = M_PI*HBAR*static_cast<double>(NX)/L;
static_assert(PMAX_HALF > 0.0);
// ***********************

/* Initial wavefunction momentum
 * NOTE: P0 > 0 makes the wave function move to the LEFT if V = 0               */
constexpr inline double P0 = -0.8*PMAX_HALF;


// ***** DON'T TOUCH *****
#define FREE_PROPAGATION 0
#define STEP             1
#define BARRIER_WELL     2
#define HARMONIC         3

// Maximum allowed value for the potential
constexpr inline double VMAX = 4.0*M_PI*HBAR/DT;
static_assert(VMAX > 0.0);
// ***********************

// Potential choice
#define POT STEP

// Central value for the potential
constexpr inline double VC = 0.5*L;

// Height of the potential barrier/well or step
constexpr inline double VH = -0.1*VMAX;

// Width of the potential barrier/well
constexpr inline double VW = 0.1*L;

// Strength of the harmonic potential
constexpr inline double K = 1.0e+05;


// Output frequency
constexpr inline size_t OUT_EVERY = 1;

// Output filename
#define FILENAME "Data.h5"

// Verbosity
constexpr inline bool VERBOSE = true;


#endif  // PARAMETERS_HH
