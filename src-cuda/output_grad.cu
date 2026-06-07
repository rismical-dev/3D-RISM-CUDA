#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_grad(double * & dulj, double * & due) {

  cout << "outputting grad to file:  " << fname + extgra << "  ..." << endl;

  ofstream out_file;
  out_file.open ((fname + extgra).c_str());

//  double dv = ce -> dv / kcal2J;
  double dv = ce -> dv;
  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    out_file << scientific 
             << setw(16) << setprecision(8) << dulj[num] * dv << " "
	     << setw(16) << setprecision(8) << dulj[num + 1] * dv << " "
	     << setw(16) << setprecision(8) << dulj[num + 2] * dv << " "
             << setw(16) << setprecision(8) << due[num] * dv << " "
	     << setw(16) << setprecision(8) << due[num + 1] * dv << " "
	     << setw(16) << setprecision(8) << due[num + 2] * dv << endl;
  }

  cout << "done." << endl;

  out_file.close();
} 
