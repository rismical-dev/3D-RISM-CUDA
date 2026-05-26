#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_qv(double * & qv) {

  cout << "outputting qv to file:  " << fname + extgra << "  ..." << endl;

  ofstream out_file;
  out_file.open ((fname + extqv).c_str());

  double dv = ce -> dv;
  int gridx = ce -> grid[0];
  int gridy = ce -> grid[1];
  int gridz = ce -> grid[2];
  for (int ig = 0; ig < ce -> ngrid; ++ig) {
     double kx = (ig % gridx - gridx / 2) * ce -> dr[0];
     double ky = ((ig / gridx) % gridy - gridy / 2) * ce -> dr[1];
     double kz = (ig / (gridx * gridy) - gridz / 2) * ce -> dr[2];
     out_file << kx << " " << ky << " " << kz << " ";
     out_file << qv[ig] * dv << endl;
  }

  cout << "done." << endl;

  out_file.close();
} 
