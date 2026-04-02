#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"


static_assert(WF == GAUSSIAN or WF == BOX);

static_assert(HBAR > 0.0);
static_assert(M    > 0.0);

static_assert(NX > 1);
static_assert(L  > 0.0);

static_assert(NT > 1);
static_assert(DT > 0.0);

static_assert(SIGMA > 0.0);
static_assert(X0 - 3.0*SIGMA > 0. and X0 + 3.0*SIGMA < L - L/static_cast<double>(NX),
              "Wave function is quite out of bounds");
static_assert(P0 >= -PMAX_HALF and P0 <= PMAX_HALF - 2.0*M_PI*HBAR/L);

static_assert(POT == FREE_PROPAGATION or POT == BARRIER_WELL or
              POT == STEP             or POT == HARMONIC);
static_assert(VC > 0.0   and VC < L);
static_assert(VH > -VMAX and VH < VMAX);
static_assert(VW > 0.0   and VW < L);
static_assert(K  > 0.0);

static_assert(OUT_EVERY > 0);
static_assert(VERBOSE or not VERBOSE);


#endif  // CHECK_PARAMETERS_HH
