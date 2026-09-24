#include <iostream>
#include <fstream>
#include <string>
#include "rism3d.h"
#include "extension.h"
#include "cuda_check.h"

void RISM3D :: write_tuv () {

  std::cout << "saving Tuv to file:  " << fname + exttuv << "  ..." << std::endl;

  // Synchronous: tuv[iv] is read on the host (written to file) right
  // below, so an async copy here would let that read race the D2H
  // transfer (same class of bug fixed in initialize_g.cu/iterate.cu).
  for (int iv = 0; iv < sv -> natv; ++iv) {
    AN_CUDA_CHECK(cudaMemcpy(tuv[iv], dt + (iv * ce -> ngrid),
			    ce -> ngrid * sizeof(double), cudaMemcpyDefault));
  }

  // Binary format (paired with read_tuv.cu): avoids both the text
  // formatting cost/precision loss of writing every double through
  // operator<< (default ostream precision is ~6 significant digits,
  // which isn't enough to restore a checkpoint bit-for-bit) and the
  // endl-per-value flush cost of the old text loop. Note: this changes
  // the .tuv file format, so a .tuv file saved by the old text-based
  // write_tuv() is no longer readable by read_tuv() -- an existing
  // checkpoint would need to be regenerated (re-run without a saved
  // file, or re-save once) after this change.
  std::ofstream out_file;
  out_file.open((fname + exttuv).c_str(), std::ios::binary);
  if (!out_file) {
    std::cerr << "Tuv file could not be opened for writing!" << std::endl;
    exit (1);
  }

  int header[4] = {ce -> grid[0], ce -> grid[1], ce -> grid[2], sv -> natv};
  out_file.write(reinterpret_cast<const char *>(header), sizeof(header));

  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file.write(reinterpret_cast<const char *>(tuv[iv]),
                    ce -> ngrid * sizeof(double));
  }

  std::cout << "done." << std::endl;

  out_file.close();
}
