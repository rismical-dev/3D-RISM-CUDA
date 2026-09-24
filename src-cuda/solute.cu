#include <cmath>
#include <limits>
#include "solute.h"
#include "cuda_check.h"

Solute :: ~Solute () {
  // Best-effort cleanup, like the other classes in this codebase (FFT3D,
  // AN2, Solvent): a destructor shouldn't abort the whole run over a
  // teardown failure, so cudaFree failures are reported rather than
  // treated as fatal. q/sig/eps/r/dq/dr are nullptr until init()/
  // setup_cuda() actually run (see the constructor), so delete[]/
  // cudaFree on any of them is always safe even if this object is
  // destroyed before those are called.
  delete[] q;
  delete[] sig;
  delete[] eps;
  delete[] r;
  cudaError_t er = cudaFree(dq);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dq) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(dr);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dr) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
}

void Solute :: init(int n) {
  num = n;
  q = new double[num];
  sig = new double[num];
  eps = new double[num];
  r = new double[num * 3];
}

// Writes the centering shift into shift[0..2], which the caller owns
// (RISM3D::read_input passes ce->shift, the array Cell's own constructor
// already allocated) instead of allocating and returning a fresh
// double[3]. The old version's caller reassigned ce->shift to this
// function's return value, silently leaking the array Cell's constructor
// had already made for it.
void Solute :: centering(double * shift) {
  double xmin, ymin, zmin, xmax, ymax, zmax;

  xmin = ymin = zmin = std::numeric_limits<double>::max();
  xmax = ymax = zmax = std::numeric_limits<double>::lowest();

  for (int n = 0; n < num; ++n) {
    int i = n * 3;
    if (xmin > r[i]) xmin = r[i];
    if (ymin > r[i + 1]) ymin = r[i + 1];
    if (zmin > r[i + 2]) zmin = r[i + 2];
    if (xmax < r[i]) xmax = r[i];
    if (ymax < r[i + 1]) ymax = r[i + 1];
    if (zmax < r[i + 2]) zmax = r[i + 2];
  }

  shift[0] = round(- (xmax - xmin) / 2 - xmin);
  shift[1] = round(- (ymax - ymin) / 2 - ymin);
  shift[2] = round(- (zmax - zmin) / 2 - zmin);

  for (int n = 0; n < num; ++n) {
    int i = n * 3;
    r[i] += shift[0];
    r[i + 1] += shift[1];
    r[i + 2] += shift[2];
  }
}

void Solute :: zero() {
  for (int n = 0; n < num; ++n) {
    q[n] = 0.0;
  }
}

void Solute :: setup_cuda() {
  AN_CUDA_CHECK(cudaMalloc(&dq, num * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&dr, num * sizeof(double3)));
  cudaMemcpyAsync(dq, q, num * sizeof(double), cudaMemcpyDefault);
  cudaMemcpyAsync(dr, r, num * sizeof(double3), cudaMemcpyDefault);
  AN_CUDA_CHECK(cudaGetLastError());
}
