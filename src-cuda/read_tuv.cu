#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"
#include "cuda_check.h"

void RISM3D :: read_tuv () {

  // Binary format, paired with write_tuv.cu's ofstream::write -- see the
  // comment there. A .tuv file saved before that change (plain text) is
  // not readable here; iterate.cu's "saved" check would find the file
  // but this would read garbage from it, so any pre-existing checkpoint
  // needs to be regenerated after this change.
  std::ifstream in_file;
  in_file.open((fname + exttuv).c_str(), std::ios::binary);
  if (!in_file) {
    std::cerr << "Tuv file could not be opened for reading!" << std::endl;
    exit (1);
  }

  int header[4];
  in_file.read(reinterpret_cast<char *>(header), sizeof(header));
  int ng0 = header[0];
  int ng1 = header[1];
  int ng2 = header[2];
  int natvo = header[3];
  if (ng0 != ce -> grid[0] || ng1 != ce -> grid[1] ||
      ng2 != ce -> grid[2] || natvo != sv -> natv) {
    std::cout << " RXRISM:  actual   NGr = " << ce -> ngrid
	 << "   NatV = " << sv -> natv << std::endl
	 << "           saved   NGr = " << ng0 * ng1 * ng2
	 << "   NatV = " << natvo << std::endl;
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    in_file.read(reinterpret_cast<char *>(tuv[iv]), ce -> ngrid * sizeof(double));
  }

  in_file.close();

  for (int iv = 0; iv < sv -> natv; ++iv) {
    AN_CUDA_CHECK(cudaMemcpyAsync(dt + (iv * ce -> ngrid), tuv[iv],
			    ce -> ngrid * sizeof(double), cudaMemcpyDefault));
  }
}
