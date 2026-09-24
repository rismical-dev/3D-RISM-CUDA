#include <iostream>
#include <fstream>
#include <unistd.h>
#include "rism3d.h"
#include "cuda_check.h"

int main (int argc, char * argv[]) {
  RISM3D * rism;
  int ch;
  int cu, dn;
  std::string input;
  std::string structure;
  std::string hs;
  std::string esp;
  bool centering = true;
  bool q0 = false;

  cu = dn = 0;
  rism = new RISM3D;

  while ((ch = getopt(argc, argv, "c:d:i:s:r:1:2:e:fz")) != -1) {
    switch (ch){
    case 'c':
      cu = atoi(optarg);
      break;
    case 'd':
      dn = atoi(optarg);
      break;
    case 'i':
      input = optarg;
      break;
    case 's':
      structure = optarg;
      break;
    case 'r':
      hs = optarg;
      break;
    case '1':
      rism -> set_ad (atof(optarg), 1);
      break;
    case '2':
      rism -> set_ad (atof(optarg), 2);
      break;
    case 'e':
      esp = optarg;
    case 'f':
      centering = false;
      break;
    case 'z':
      q0 = true;
      break;
    }
  }

  if (input.empty() || structure.empty()) {
    if (argv[optind] == NULL) {
      std::cout << "No input file!" << std::endl;
      return (1);
    }
    input = argv[optind];
  }

  std::cout << "Set device " << dn << std::endl;
  AN_CUDA_CHECK(cudaSetDevice(dn));
  if (cu > 0) std::cout << "Charge up " << cu << std::endl;
  rism -> initialize(input, structure, esp, hs, centering, q0);
  rism -> iterate(cu);
  rism -> output();
  delete rism;

  return(0);
}
