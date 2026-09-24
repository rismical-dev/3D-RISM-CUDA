#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

// xmu/xmu2/euv are only read here, never reassigned, so they're taken as
// plain const double* rather than double * & (a reference to the caller's
// pointer, which this function has no reason to need).
void RISM3D :: output_xmu(const double * xmu, const double * xmu2,
                          const double * euv,
                          double dft, double pmv, double pressure) {

  std::string fxmu;
  if (!zero) {
    fxmu = fname + extxmu;
  } else {
    fxmu = fname + extxmu0;
  }
  std::ofstream out_file;
  out_file.open(fxmu.c_str());
  if (!out_file) {
    std::cerr << "Xmu file could not be opened!" << std::endl;
    exit(1);
  }

  double ibeta = avogadoro * boltzmann * sv -> temper;

  double pcterm = - pressure * pmv * ibeta;

  out_file << " $RESULT" << std::endl;

  if (rmdft) {
    out_file << "SFE_RMDFT= " << std::fixed << std::setprecision(5)
  	   << ibeta * dft << " !(J/mol)" << std::endl;
    out_file << std::endl;
  }

  double xmua = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    xmua += xmu[iv];
  }

  out_file << "SFE_SC= " << std::fixed << std::setprecision(5)
  	   << ibeta * xmua << " !(J/mol)" << std::endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  SFEC_SC(" << iv << ")= " << std::fixed << std::setprecision(5)
             << ibeta * xmu[iv] << std::endl;
  }
  out_file << std::endl;

  if (clos == 0) {
    xmua = 0.0;
    for (int iv = 0; iv < sv -> natv; ++iv) {
      xmua += xmu2[iv];
    }

    out_file << "SFE_SC_HNC= " << std::fixed << std::setprecision(5)
    	     << ibeta * xmua << " !(J/mol)" << std::endl;

    for (int iv = 0; iv < sv -> natv; ++iv) {
      out_file << "  SFEC_SC_HNC(" << iv << ")= " << std::fixed << std::setprecision(5)
               << ibeta * xmu2[iv] << std::endl;
    }
    out_file << std::endl;
  }

  xmua = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    xmua += xmu[sv -> natv + iv];
  }

  out_file << "SFE_GF= " << std::fixed << std::setprecision(5)
	   << ibeta * xmua << " !(J/mol)" << std::endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  SFEC_GF(" << iv << ")= " << std::fixed << std::setprecision(5)
             << ibeta * xmu[sv -> natv + iv] << std::endl;
  }
  out_file << std::endl;

  double dv = ce -> dv;
  xmua = 0.0;
  for (int i = 0; i < su -> num * sv -> natv * 2; ++i) {
    xmua += euv[i];
  }
  out_file << "SE= " << std::fixed << std::setprecision(5)
	   << xmua * dv << " !(J/mol)" << std::endl;

  xmua = 0.0;
  for (int i = 0; i < su -> num * sv -> natv; ++i) {

    xmua += euv[i];
  }
  out_file << "  SE_LJ= " << std::fixed << std::setprecision(5)
           << xmua * dv << std::endl;

  xmua = 0.0;
  for (int i = su -> num * sv -> natv; i < su -> num * sv -> natv * 2; ++i) {
    xmua += euv[i];
  }
  out_file << "  SE_ES= " << std::fixed << std::setprecision(5)
           << xmua * dv << std::endl;
  out_file << std::endl;

  out_file << "PMV= " << std::fixed << std::setprecision(5)
           << pmv * 1000 << " !(L/mol)" << std::endl;

  out_file << "Pressure= " << std::fixed << std::setprecision(5)
           << ibeta * pressure << " !(J/m^3)" << std::endl;

  out_file << "Correction_Term= " << std::fixed << std::setprecision(5)
             << pcterm << " !(J/mol)" << std::endl;

  out_file << " $END" << std::endl;

  out_file.close();
}
