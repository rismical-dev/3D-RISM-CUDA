#ifndef EXCP_FD_H
#define EXCP_FD_H

// Scaled-particle-theory excess chemical potential correction for a hard
// sphere of packing fraction pfhs in a fluid of number density dens.
// Defined in excp_fd.cc; declared here so callers (currently cal_rmdft.cu)
// don't each carry their own local forward declaration that could drift
// out of sync with the definition's signature.
double excp_fd (double pfhs, double dens);

#endif  // EXCP_FD_H
