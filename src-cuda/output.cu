#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);

  if (outlist.find("m") != string::npos) {
    double pmv = cal_pmv();
    double pressure = cal_pressure();
    double * xmu = new double[sv -> natv * 2];
    double * xmu2 = new double[sv -> natv];
    double dft;

    cal_exchem(xmu, xmu2);
    if (rmdft) dft = cal_rmdft();
    output_xmu(xmu, xmu2, dft, pmv, pressure);
    delete[] xmu, xmu2;
  }

  if (outlist.find("d") != string::npos) {
    double * du;
    du = new double[su -> num * 3];
    cal_grad(du);
    output_grad(du);
    delete[] du;
  }

  if (outlist.find("c") != string::npos) {
    output_cuv();
  }

  if (outlist.find("g") != string::npos) {
    output_guv();
  }

  if (outlist.find("h") != string::npos) {
    output_huv();
  }

  if (outlist.find("e") != string::npos) {
    double euv = cal_euv();
    output_euv(euv);
  }
} 
