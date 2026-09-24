#ifndef SPLINE_H
#define SPLINE_H
#include <vector>

// Natural cubic spline fit/evaluation, shared by Solvent::spline() and
// Solvent::spline2() (solvent_spline.cu / solvent_spline2.cu). Declared
// here instead of being locally forward-declared in each caller, so the
// two files can't drift out of sync with spline.cc/splint.cc's actual
// signatures.
//
// xp/yp are read-only (never reassigned), so they're passed as plain
// const double* rather than the original double*& -- that reference-to-
// pointer added indirection without ever being used to hand a different
// pointer back to the caller. coe is likewise read-only in splint() and
// only has its pointees written through (never reassigned or resized)
// in spline(), so both take it as a const reference.
void spline (const double * xp, const double * yp, int np,
             const std::vector <double *> & coe);
double splint (const double * xp, const double * yp,
               const std::vector <double *> & coe, int np, double x);

#endif  // SPLINE_H
