#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_huv() {

  std::cout << "outputting Huv to file:  " << fname + exthuv << "  ..." << std::endl;

  std::ofstream out_file;
  out_file.open((fname + exthuv).c_str());
  if (!out_file) {
    std::cerr << "Huv file could not be opened!" << std::endl;
    exit(1);
  }

  out_file << ce -> box[0] << " "
	   << ce -> box[1] << " "
	   << ce -> box[2] << " "
	   << std::endl
	   << ce -> grid[0] << " "
	   << ce -> grid[1] << " "
	   << ce -> grid[2] << " "
	   << std::endl
	   << su -> num << std::endl;

  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    out_file << su -> q[iu] << " "
	     << su -> sig[iu] << " "
	     << su -> eps[iu] << " "
	     << su -> r[num] << " "
	     << su -> r[num + 1] << " "
	     << su -> r[num + 2] << std::endl;
  }
  out_file << sv -> natv << std::endl;

  for (int ig = 0; ig < ce -> ngrid; ++ig) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      out_file << huv[iv][ig].real() << " ";
    }
    out_file << std::endl;
  }

  std::cout << "done." << std::endl;

  out_file.close();
}
