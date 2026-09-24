#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_guv() {

  std::cout << "outputting Guv to file:  " << fname + extguv << "  ..." << std::endl;

  std::ofstream out_file((fname + extguv).c_str(), std::ios::out | std::ios::binary);
  if (!out_file) {
    std::cerr << "Guv file could not be opened!" << std::endl;
    exit(1);
  }
  out_file.write((char *) &(sv -> natv), sizeof(int));
  out_file.write((char *) &(ce -> grid[0]), sizeof(int) * 3);
  out_file.write((char *) &(ce -> box[0]), sizeof(double) * 3);
  out_file.write((char *) &(ce -> shift[0]), sizeof(double) * 3);
  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      double tmp = guv[iv][ig].real();
      out_file.write((char *) &tmp, sizeof(double));
    }
  }
  out_file.close();

  std::cout << "done." << std::endl;
}
