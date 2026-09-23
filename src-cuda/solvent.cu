#include "solvent.h"
#include "cuda_check.h"

// Solvent::~Solvent() was declared in solvent.h but never defined anywhere
// in this codebase. That went unnoticed because RISM3D's own destructor
// had a comma-operator bug (see rism3d.h) that meant "delete sv;" never
// actually ran -- so the linker never needed to resolve this symbol.
// Fixing that bug surfaced the missing definition. It needs its own .cu
// file (rather than living in solvent_read.cc, a plain .cc compiled with
// g++ and no CUDA runtime headers on its include path, or being folded
// into one of the spline .cu files that also happen to have other jobs)
// because it has to cudaFree the device buffers this class owns.
Solvent :: ~Solvent () {
  void dealloc2D (vector <double *> &);
  void dealloc3D (vector <vector <double *> > &);

  // ttab / ttab2 are NOT freed here. They are already delete[]'d
  // unconditionally, inline, at the end of spline() / spline2() once
  // those functions finish consuming them (see solvent_spline.cu /
  // solvent_spline2.cu). Nothing marks a raw pointer as "already freed",
  // so deleting them again here would double-free -- the same hazard
  // AN2 and FFT3D originally had for their own buffers, except there the
  // fix was "the owner deletes it exactly once, elsewhere, on purpose"
  // rather than "never actually allocated".
  delete[] sigv;
  delete[] epsv;
  delete[] qv;
  delete[] rhov;
  delete[] pfhs;
  delete[] wfk0;

  // xvva/cvva/wfka are the host-side staging buffers spline()/spline2()
  // build before copying them to the device (dx/dc/dw below); nothing
  // ever freed the host copies once that copy was made. xvv/cvv/chs/wfk
  // are the raw tabulated-data buffers read() allocates: xvv is already
  // emptied by spline()'s own dealloc3D(xvv), so this is a safe no-op
  // for it, but cvv/chs/wfk (only ever populated when a hard-sphere file
  // is given) were never freed anywhere -- a genuine, separate leak from
  // the missing destructor itself. dealloc2D/dealloc3D are safe to call
  // on an empty vector (nothing was ever allocated, or it already was
  // deallocated), so none of these need an rmdft/hs-string guard here.
  dealloc3D(xvva);
  dealloc3D(cvva);
  dealloc3D(xvv);
  dealloc3D(cvv);
  dealloc2D(wfka);
  dealloc2D(chs);
  dealloc2D(wfk);

  // dx is always allocated (spline() runs unconditionally); dc/dw/drho
  // only when rmdft is true (spline2() runs). When it isn't, they stay
  // nullptr from the constructor, and cudaFree(nullptr) is a documented
  // no-op, so no extra guard is needed here either. Best-effort cleanup,
  // like FFT3D's and AN2's destructors: report a failure instead of
  // aborting the whole run over it.
  cudaError_t er = cudaFree(dx);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dx) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(dc);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dc) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(dw);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dw) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(drho);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(drho) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
}
