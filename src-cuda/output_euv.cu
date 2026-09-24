#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

// euv is only read here, never reassigned, so it's taken as plain
// const double* rather than double * & (a reference to the caller's
// pointer, which this function has no reason to need).
void RISM3D :: output_euv (const double * euv) {

  std::cout << "outputting euv to file:  " << fname + exteuv << "  ..." << std::endl;

  std::ofstream out_file;
  out_file.open ((fname + exteuv).c_str());
  if (!out_file) {
    std::cerr << "Euv file could not be opened!" << std::endl;
    exit(1);
  }

  double dv = ce -> dv;
  // The ES half of euv[] starts right after the LJ half (su->num * sv->natv
  // entries) -- matches cal_euv.cu's i=1 term in
  // e[(natv*num)*i + num*iv+iu]. Named for what it is, rather than reusing
  // "i" as in that 0/1 block index, and hoisted out of the iv loop since it
  // doesn't depend on iv.
  const int block_offset = su -> num * sv -> natv;
  for (int iu = 0; iu < su -> num; ++iu) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      out_file << std::scientific << std::setw(16) << std::setprecision(8)
               << euv[su -> num * iv + iu] * dv << " ";
    }
    for (int iv = 0; iv < sv -> natv; ++iv) {
      out_file << std::scientific << std::setw(16) << std::setprecision(8)
               << euv[su -> num * iv + block_offset + iu] * dv << " ";
    }
    out_file << std::endl;
  }

  std::cout << "done." << std::endl;

  out_file.close () ;
}
