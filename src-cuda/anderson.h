#ifndef AN2_H
#define AN2_H

#include "cell.h"
#include "solvent.h"

class AN2 {
public:
  AN2 () : m(0.0), mp(0.0), count(0), dtp(nullptr), drp(nullptr),
           irho(nullptr), s(nullptr), ds(nullptr), a(nullptr), c(nullptr),
           x(nullptr), s1(0.0), s2(0.0), niv(0), ngrid(0), binary(0),
           initialized(false) {}
  ~AN2 ();

  // AN2 owns device buffers (dtp, drp, ds) and host buffers (s, a, c, x,
  // irho), all allocated in initialize(). Copying it would give two
  // objects the same pointers, and whichever is destroyed first would
  // free memory the other still thinks it owns (the second destructor
  // call would then double-free). Nothing in the codebase copies an AN2
  // today -- RISM3D only ever holds one through a pointer -- so disabling
  // copy costs nothing and rules that bug class out entirely.
  AN2 (const AN2 &) = delete;
  AN2 & operator= (const AN2 &) = delete;

  void initialize (Cell *, Solvent *);
  void calculate (double *, double *);
  double m;
  double mp;
  int count;
private:
  void cal_theta (double *, double *);
  double * dtp;
  double * drp;
  double * irho;
  double * s;
  double * ds;
  double * a;
  double * c;
  double * x;
  double s1;
  double s2;
  int niv;
  int ngrid;
  int binary;
  dim3 g, b;
  bool initialized;
};

#endif // AN2_H
