#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"

static_assert(WF == GAUSSIAN or WF == BOX);
static_assert(HBAR > 0.0);
static_assert(NX > 1);
static_assert(L > 0.0);
static_assert(SIGMA > 0.0);
static_assert(X0 - 3.0*SIGMA > -L/2.0 and X0 + 3.0*SIGMA < L/2.0 - L/static_cast<double>(NX),
              "Wave function is quite out of bounds");
static_assert(P0 >= -PMAX_HALF and P0 <= PMAX_HALF - 2.0*M_PI*HBAR/L);

#endif  // CHECK_PARAMETERS_HH
