#ifndef SOLUTE_H
#define SOLUTE_H

class Solute {
 public:
  Solute () : q(nullptr), sig(nullptr), eps(nullptr), r(nullptr),
              dq(nullptr), dr(nullptr), num(0) {}
  ~Solute ();

  // Solute owns host arrays (q, sig, eps, r) and device buffers (dq, dr).
  // Copying it would give two objects the same pointers, and whichever
  // is destroyed first would free memory the other still thinks it
  // owns. Nothing in the codebase copies a Solute today (RISM3D only
  // ever holds one through a pointer), so disabling copy costs nothing
  // and rules that bug class out entirely.
  Solute (const Solute &) = delete;
  Solute & operator= (const Solute &) = delete;

  void init (int);
  void centering (double *);
  void zero ();
  void setup_cuda ();
  double * q;
  double * sig;
  double * eps;
  double * r;
  double * dq;
  double3 * dr;
  int num;
};

#endif
