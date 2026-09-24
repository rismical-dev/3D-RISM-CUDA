#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

// qv is only read here, never reassigned, so it's taken as plain
// const double* rather than double * & (a reference to the caller's
// pointer, which this function has no reason to need).
void RISM3D :: output_qv(const double * qv) {

  std::cout << "outputting qv to file:  " << fname + extqv << "  ..." << std::endl;

  std::ofstream out_file;
  out_file.open ((fname + extqv).c_str());
  if (!out_file) {
    std::cerr << "Qv file could not be opened!" << std::endl;
    exit(1);
  }

  double dv = ce -> dv;
  int gridx = ce -> grid[0];
  int gridy = ce -> grid[1];
  int gridz = ce -> grid[2];
  for (int ig = 0; ig < ce -> ngrid; ++ig) {
    double kx = (ig % gridx - gridx / 2) * ce -> dr[0];
    double ky = ((ig / gridx) % gridy - gridy / 2) * ce -> dr[1];
    double kz = (ig / (gridx * gridy) - gridz / 2) * ce -> dr[2];
    out_file << "    "
             << std::fixed
             << std::setw(12) << std::setprecision(4) << kx
             << std::setw(12) << std::setprecision(4) << ky
             << std::setw(12) << std::setprecision(4) << kz
             << "  "
             << std::scientific
             << std::setw(16) << std::setprecision(8) << qv[ig] * dv
             << std::endl;
  }

  std::cout << "done." << std::endl;

  out_file.close();
}
