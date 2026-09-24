#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

// ad is only read here, never reassigned, so it's taken as plain
// const double* rather than double * & (a reference to the caller's
// pointer, which this function has no reason to need).
void RISM3D :: output_ad(const double * ad) {

  std::stringstream ss;
  ss << fname << "-" << adswitch << "-" << lambda;

  std::cout << "outputting ad to file:  " << ss.str() << "  ..." << std::endl;

  std::ofstream out_file;
  out_file.open (ss.str().c_str());
  if (!out_file) {
    std::cerr << "Ad file could not be opened!" << std::endl;
    exit(1);
  }

  double dv = ce -> dv;
  for (int iu = 0; iu < su -> num; ++iu) {
    out_file << std::scientific
             << std::setw(16) << std::setprecision(8) << ad[iu] * dv << std::endl;
  }

  std::cout << "done." << std::endl;

  out_file.close();
}
