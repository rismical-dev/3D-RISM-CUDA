#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_solvent (string hs) {
  if (hs != "") rmdft = true;
  sv -> read(fsolvent, hs);
  sv -> spline(ga, indga, nga, ce -> ngrid, rmdft);
  if (rmdft) sv -> spline2(ga, indga, nga, ce -> ngrid);
}
