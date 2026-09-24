#include <iostream>
#include "rism3d.h"

// calculate and return P/kBT as pressure

double RISM3D :: cal_pressure () {

  // sv -> rhov[] is in molecules/nm^3 (see Solvent::read(), solvent_read.cc:
  // "rhov[i] *= avogadoro * 1.0e-27" converts a molar concentration to
  // that unit); this converts back up to molecules/m^3 for the P/kBT
  // result documented below.
  constexpr double NM3_TO_M3 = 1.0e30;

  double pressure = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    pressure +=  sv -> rhov[iv];
  }

  pressure *= NM3_TO_M3 / avogadoro;

  double ibeta = avogadoro * boltzmann * sv -> temper;
  pressure = 0.5 * (pressure + 1.0 / (sv -> xt * ibeta));
// P/kBT (pressure ) is in [mol/m^3]

  return pressure;
}
