#ifndef SOLVENT_H
#define SOLVENT_H
#include <stdlib.h>
#include <vector>
#include <string>
using namespace std;

class Solvent {
 public:
  Solvent () : dx(nullptr), dc(nullptr), dw(nullptr), drho(nullptr),
               rhov(nullptr), qv(nullptr), sigv(nullptr), epsv(nullptr),
               pfhs(nullptr), wfk0(nullptr), temper(0.0), xt(0.0), natv(0),
               ttab(nullptr), ttab2(nullptr), ntab(0), ntab2(0) {}
  ~Solvent ();
  void read (string, string);
  void spline (vector <double> &, int * &, int, int, bool);
  void spline2 (vector <double> &, int * &, int, int);
  vector <vector <double *> > xvva;
  vector <vector <double *> > cvva;
  vector <double *> wfka;
  double * dx;
  double * dc;
  double * dw;
  double * drho;
  double * rhov;
  double * qv;
  double * sigv;
  double * epsv;
  double * pfhs;
  double * wfk0;
  double temper;
  double xt;
  int natv;
 private:
  vector <vector <double *> > xvv;
  vector <vector <double *> > cvv;
  vector <double *> chs;
  vector <double *> wfk;
  double * ttab;
  double * ttab2;
  int ntab;
  int ntab2;
};

#endif
