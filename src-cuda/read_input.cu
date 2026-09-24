#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "version.h"

bool isNaturalNumber(const std::string& str) {
  try {
    size_t pos;
    int num = std::stoi(str, &pos);
    if (pos != str.length()) {
      return false;
    }
    return num >= 0;
  } catch (const std::exception& e) {
    return false;
  }
}

void RISM3D :: read_input (std::string control, std::string structure, bool centering) {

  std::ifstream in_file;
  in_file.open (control.c_str());
  if (!in_file) {
    std::cerr << "Control file could not be opened!" << std::endl;
    exit (1);
  }

  std::cout << "reading input data file  : " << control << std::endl;

  std::string check;
  in_file >> outlist >> co -> ksave >> check;
  if (check != version) {
    std::cout << "Error: This input file is for old version." << std::endl;
    exit (1);
  }

  std::string closure;
  in_file >> closure;
  if (closure == "KH") {
    clos = 0;
  } else if (closure == "HNC") {
    clos = 1;
  } else {
    std::cout << "Error: Unexpected closure switch " << std::endl;
    exit(1);
  }
  in_file >> fsolvent;
  in_file >> co -> convergence >> co -> maxstep;
  in_file >> ma -> count >> ma -> m >> ma -> mp;
  in_file >> ce -> box[0] >> ce -> box[1] >> ce -> box[2];
  in_file >> ce -> grid[0] >> ce -> grid[1] >> ce -> grid[2];

  if (!structure.empty()) {
    in_file.close ();
    in_file.open (structure.c_str());
    if (!in_file) {
      std::cerr << "Structure file could not be opened!" << std::endl;
      exit (1);
    }
    std::cout << "reading solute data file : " << structure << std::endl;
  }

  int num;

  std::string tmp;
  in_file >> tmp;

  if (isNaturalNumber(tmp)) {
    num = std::stoi(tmp);
  } else {
    std::cout << "Error: Number of solute atom is not a natural number." << std::endl;
    exit(1);
  }

  ce -> setup();
  su -> init(num);

  for (int iu = 0; iu < su -> num; ++iu) {
    int n = iu * 3;
    in_file >> su -> sig[iu] >> su -> eps[iu] >> su -> q[iu]
	    >> su -> r[n] >> su -> r[n + 1] >> su -> r[n + 2];
  }

  in_file.close ();

  if (centering) {
    su -> centering(ce -> shift);
  }

  if (zero) {
    su -> zero();
  }

  su -> setup_cuda();
}
