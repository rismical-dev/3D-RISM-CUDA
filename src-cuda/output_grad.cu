#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

// dulj/due are only read here, never reassigned, so they're taken as plain
// const double* rather than double * & (a reference to the caller's
// pointer, which this function has no reason to need).
void RISM3D :: output_grad(const double * dulj, const double * due) {

  std::cout << "outputting grad to file:  " << fname + extgra << "  ..." << std::endl;

  std::ofstream out_file;
  out_file.open ((fname + extgra).c_str());
  if (!out_file) {
    std::cerr << "Grad file could not be opened!" << std::endl;
    exit(1);
  }

  double dv = ce -> dv;
  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    out_file << std::scientific
             << std::setw(16) << std::setprecision(8) << dulj[num] * dv << " "
	     << std::setw(16) << std::setprecision(8) << dulj[num + 1] * dv << " "
	     << std::setw(16) << std::setprecision(8) << dulj[num + 2] * dv << " "
             << std::setw(16) << std::setprecision(8) << due[num] * dv << " "
	     << std::setw(16) << std::setprecision(8) << due[num + 1] * dv << " "
	     << std::setw(16) << std::setprecision(8) << due[num + 2] * dv << std::endl;
  }

  std::cout << "done." << std::endl;

  out_file.close();
}
