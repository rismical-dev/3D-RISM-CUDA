#ifndef FFT3D_H
#define FFT3D_H
#include <cufft.h>
#include "cell.h"

class FFT3D {
 public:
  // Direction codes for execute()'s "key" argument, replacing the bare
  // -1 / 1 magic numbers that used to appear at every call site.
  // Values are kept numerically identical to the old convention (and to
  // plain int) so existing callers that still pass literal -1 / 1
  // (calculate.cu, cal_rmdft.cu) keep compiling unchanged; they can be
  // migrated to FFT3D::FORWARD / FFT3D::BACKWARD incrementally.
  static const int FORWARD  = -1;
  static const int BACKWARD =  1;

  FFT3D () : plan(0), dkf(nullptr), dir(nullptr),
             volf(0.0), volb(0.0), ngrid(0), initialized(false) {}
  ~FFT3D ();

  // FFT3D owns a cuFFT plan and two device buffers. Copying it would give
  // two objects the same plan handle and the same device pointers, so
  // whichever one is destroyed first would free resources the other still
  // thinks it owns (and the second destructor call would double-free).
  // Nothing in the codebase copies an FFT3D today -- RISM3D only ever
  // holds one through a pointer -- so disabling copy costs nothing and
  // just rules that bug out for good. Move is left disabled too, for the
  // same reason: there is no code that needs it, and adding it correctly
  // only matters once something does.
  FFT3D (const FFT3D &) = delete;
  FFT3D & operator= (const FFT3D &) = delete;

  void initialize (Cell *);
  void execute (double2 *, int key);
 private:
  cufftHandle plan;
  dim3 g, b;
  double2 * dkf;
  int * dir;
  double volf, volb;
  int ngrid;
  bool initialized;
};
#endif
