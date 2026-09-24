#ifndef RISM3D_H
#define RISM3D_H

#include <complex>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "physical.h"
#include "cell.h"
#include "control.h"
#include "solute.h"
#include "solvent.h"
#include "anderson.h"
#include "fft3d.h"

class RISM3D {
public:
  RISM3D () {ce = new Cell; co = new Control; su = new Solute;
    sv = new Solvent; ma = new AN2; fft = new FFT3D;}
  // ma is deliberately NOT deleted here: iterate() already does
  // "delete ma;" itself once the SCF loop converges (AN2 isn't needed
  // again after that), so deleting it a second time here would be a
  // double-free. fft has no such early delete -- iterate() keeps it
  // alive on purpose (cal_rmdft(), called from output(), still needs
  // fft->execute()) -- and previously wasn't deleted anywhere at all,
  // which was a real leak; it belongs here, once output() has returned.
  ~RISM3D () {delete ce; delete co; delete su; delete sv; delete fft;}

  // RISM3D owns several heap-allocated sub-objects (ce, co, su, sv, ma,
  // fft), some of which (ma, fft) in turn own CUDA/cuFFT resources.
  // Copying a RISM3D would give two objects the same pointers, and
  // whichever is destroyed first would free state the other still thinks
  // it owns. Nothing in the codebase copies a RISM3D today (main.cu only
  // ever holds one through a pointer), so disabling copy costs nothing
  // and rules that bug class out entirely.
  RISM3D (const RISM3D &) = delete;
  RISM3D & operator= (const RISM3D &) = delete;

  void set_ad (double, int);
  void initialize (std::string, std::string, std::string, std::string, bool, bool);
  void iterate (int);
  void output ();
private:
  void cal_ad1 (double * &);
  void cal_ad2 (double * &);
  void add_tuv (double);
  void cal_Coulomb (std::string);
  void cal_euv (double * &);
  void cal_exchem (double * &, double * &);
  void cal_qv (double * &);
  double cal_rmdft ();
  void cal_grad (double * &, double * &);
  void cal_LJ ();
  double cal_pmv ();
  void cal_potential (std::string);
  double cal_pressure ();
  double cal_rms ();
  void calculate (double);
  void initialize_g ();
  void initialize_tuv ();
  void output_ad (const double *);
  void output_cuv ();
  void output_euv (const double *);
  void output_grad (const double *, const double *);
  void output_guv ();
  void output_huv ();
  void output_qv (const double *);
  void output_xmu (const double *, const double *, const double *, double, double, double);
  void read_input (std::string, std::string, bool);
  void read_tuv ();
  void set_fname (std::string, std::string);
  void set_cuda ();
  void set_solvent (std::string);
  void write_tuv ();

  // Host
  std::vector <std::complex <double> *> guv;
  std::vector <std::complex <double> *> huv;
  std::vector <double *> tuv;
  std::vector <double> ga;
  // siguv/epsuv used to live here as members, but they were only ever
  // used inside cal_LJ() as staging buffers for the host-side sig/eps
  // combination rules before copying them to the device -- nothing else
  // in the codebase touched them, and nothing ever freed them. They're
  // now plain locals in cal_LJ.cu, deleted right after use.
  double lambda;
  int * indga;
  int clos;
  int nga;
  int adswitch = 0;
  std::string outlist;
  std::string fsolvent;
  std::string fname;
  dim3 g, b;
  bool rmdft = false;
  bool zero = false;

  // Device
  double2 * dguv;
  double2 * dhuv;
  double * dt;
  double * dtr;
  double * du;
  double * de;
  double3 * dgv;
  double * dsig;
  double * deps;
  double * dfr;
  double2 * dfk;
  double * ds;
  int * dindga;

  Cell * ce;
  Control * co;
  Solute * su;
  Solvent * sv;
  AN2 * ma;
  FFT3D * fft;
};

#endif
