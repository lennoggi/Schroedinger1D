#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"

static_assert(WF == GAUSSIAN or WF == BOX);
static_assert(HBAR > 0.);
static_assert(N > 1);
static_assert(L > 0.);
static_assert(SIGMA > 0.);
static_assert(X0 - 3.*SIGMA > -L/2. and X0 + 3*SIGMA < L/2. - L/static_cast<double>(N),
              "Wave function is quite out of bounds");
static_assert(P0 >= -PMAX_HALF and P0 <= PMAX_HALF - 2.*M_PI*HBAR/(1.*L));

#endif  // CHECK_PARAMETERS_HH
