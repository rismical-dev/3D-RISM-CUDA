#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "solvent.h"
#include "physical.h"

void Solvent :: read(string fsolvent, string hs) {
  void alloc2D (vector <double *> &, int, int);
  void alloc3D (vector <vector <double * > > &, int, int, int);

  ifstream in_file;
  in_file.open(fsolvent.c_str());
  if (!in_file) {
    std::cerr << "Solvent file could not be opened!" << std::endl;
    exit (1);
  }

  double * sig;
  double * eps;
  double * q;
  double * den;
  int * sol;
  string dummy;  
  double dr;
  int num;
  auto i = 0;
  std::string line;
  while (std::getline(in_file, line)) {
    std::istringstream iss(line);
    if (i == 1) {
      iss  >> dummy >> num >> natv >> ntab >> dr >> dummy;
      alloc3D (xvv, natv, natv, ntab);
      sig = new double[num];
      eps = new double[num];
      q = new double[num];
      den = new double[num];
      sol = new int[num];
    }
    if (i == 2) {
      iss  >> dummy >> dummy >> temper >> xt;
    }
    if (i > 3 && i <= num + 3) {
      int i2 = i - 4;
      iss >> dummy >> dummy >> dummy >> sol[i2]
          >> sig[i2] >> eps[i2] >> q[i2]
          >> dummy >> dummy >> dummy >> den[i2];
    }
    if (i > num + 3 && i < num + 4 + ntab) {
      int i2 = i - num - 4;
      for (int iv2 = 0; iv2 < natv; ++iv2) {
	for (int iv1 = 0; iv1 < natv; ++iv1) {
	  in_file >> xvv[iv2][iv1][i2];
	}
      }
    }
    ++i;
  }
  in_file.close();

  sigv = new double[natv];
  epsv = new double[natv];
  qv = new double[natv];
  rhov = new double[natv];

  for (int i = 0; i < natv; ++i) {
    rhov[i] = 0.0;
  }
  for (int i = 0; i < num; ++i) {
    int i2 = abs(sol[i]) - 1;
    sigv[i2] = sig[i];
    epsv[i2] = eps[i];
    qv[i2] = q[i];
    rhov[i2] += den[i];
  }
  for (int i = 0; i < natv; ++i) {
    rhov[i] *= (avogadoro * 1.0e-27);
  }

  double dk = M_PI / (dr * ntab);
  ttab = new double[ntab];
  for (int i = 0; i < ntab; ++i) {
    ttab[i] = dk * (i + 1);
  }

  if (hs != "") {
    pfhs = new double[natv]{};
    wfk0 = new double[natv]{};
    alloc3D (cvv, natv, natv, ntab);

    size_t lastdp = fsolvent.rfind('.');
    std::string cvk = fsolvent.substr(0, lastdp) + ".cvk";

    in_file.open(cvk.c_str());
    if (!in_file) {
      std::cerr << "cvk file could not be opened!" << std::endl;
      exit (1);
    }

    i = 0;
    std::string line;
    while (std::getline(in_file, line)) {
      std::istringstream iss(line);
      if (i > num + 3 && i < num + 4 + ntab) {
	int i2 = i - num - 4;
	for (int iv2 = 0; iv2 < natv; ++iv2) {
	  for (int iv1 = 0; iv1 < natv; ++iv1) {
	    in_file >> cvv[iv2][iv1][i2];
	  }
	}
      }
      ++i;
    }
    in_file.close();

    in_file.open(hs.c_str());
    if (!in_file) {
      std::cerr << "Hardsphare file could not be opened!" << std::endl;
      exit (1);
    }

    i = 0;
    while (std::getline(in_file, line)) {
      std::istringstream iss(line);
      if (i == 0) {
	iss >> ntab2;
	ttab2 = new double[ntab2];
	alloc2D (chs, natv, ntab2);
	alloc2D (wfk, natv, ntab2);
      }
      if (i == 1) {
	iss  >> pfhs[0] >> wfk0[0];
      }
      if (i >= 2) {
	auto i2 = i - 2;
	iss >> ttab2[i2] >> chs[0][i2] >> wfk[0][i2];
      }
      ++i;
    }
    in_file.close();
  }

  delete[] sig;
  delete[] eps;
  delete[] q;
  delete[] den;
  delete[] sol;
}
