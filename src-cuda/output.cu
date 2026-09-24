#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  std::transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);

  double * euv;
  if (outlist.find("m") != std::string::npos || outlist.find("e") != std::string::npos) {
    euv = new double[su -> num * sv -> natv * 2];
    cal_euv(euv);
  }

  if (outlist.find("m") != std::string::npos) {
    double pmv = cal_pmv();
    double pressure = cal_pressure();
    double * xmu = new double[sv -> natv * 2];
    double * xmu2 = new double[sv -> natv];
    double dft = 0.0;

    cal_exchem(xmu, xmu2);
    if (rmdft) dft = cal_rmdft();
    output_xmu(xmu, xmu2, euv, dft, pmv, pressure);
    delete[] xmu;
    delete[] xmu2;
  }

  if (outlist.find("d") != std::string::npos) {
    double * dulj;
    double * due;
    dulj = new double[su -> num * 3];
    due = new double[su -> num * 3];
    cal_grad(dulj, due);
    output_grad(dulj, due);
    delete[] dulj;
    delete[] due;
  }

  if (outlist.find("c") != std::string::npos) {
    output_cuv();
  }

  if (outlist.find("g") != std::string::npos) {
    output_guv();
  }

  if (outlist.find("h") != std::string::npos) {
    output_huv();
  }

  if (outlist.find("a") != std::string::npos) {
    double * ad;
    ad = new double[su -> num];
    if (adswitch == 1) {
      cal_ad1(ad);
    } else {
      cal_ad2(ad);
    }
    output_ad(ad);
    delete[] ad;
  }

  if (outlist.find("e") != std::string::npos) {
    output_euv(euv);
  }

  if (outlist.find("q") != std::string::npos) {
    double * qv = new double[ce -> ngrid];
    cal_qv(qv);
    output_qv(qv);
    delete[] qv;
  }
}
