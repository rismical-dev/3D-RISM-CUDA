#include <iostream>
#include <thrust/device_vector.h>
#include "anderson.h"
#include "cuda_check.h"

// Internal to this file: a hand-rolled Gauss-elimination solve for the
// small (2x2, in this class's ordinary-Anderson case) normal-equation
// system in AN2::calculate() (defined at the bottom of this file). Marked
// static (internal linkage) so its short, generic name can't collide with
// an unrelated "leg" defined in another translation unit at link time.
// (A local, block-scope forward declaration -- the trick used below for
// the __global__ kernels -- can't carry "static": C++ does not allow a
// static function to be declared inside another function, so this needs
// a real file-scope declaration instead.)
static void leg(double * a, double * x, double * b, int n, int nn);

AN2 :: ~AN2 () {
  if (!initialized) return;
  // Best-effort cleanup, like FFT3D's destructor: a destructor shouldn't
  // abort the whole run over a teardown failure, so failures are reported
  // rather than treated as fatal.
  delete[] a;
  delete[] c;
  delete[] x;
  delete[] s;
  delete[] irho;
  cudaError_t er = cudaFree(dtp);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dtp) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(drp);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(drp) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(ds);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(ds) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
}


void AN2 :: initialize (Cell * ce, Solvent * sv) {
  ngrid = ce -> ngrid;
  niv = sv -> natv;
  binary = 0;
  b.x = ce -> grid[0];
  g.x = ce -> grid[1];
  g.y = ce -> grid[2];

  // Every kernel below is launched as <<< g, b >>>, i.e. with b.x ==
  // ce->grid[0] threads per block. This used to be checked here against
  // the device's maxThreadsPerBlock, but that check is unreachable now:
  // RISM3D::set_cuda() (set_cuda.cu) checks the exact same ce->grid[0]
  // against the exact same device limit, and always runs first --
  // set_cuda() is called from RISM3D::initialize(), while this function is
  // only reached later, from iterate()'s "ma -> initialize(ce, sv)". If
  // grid[0] were too large, set_cuda() would already have exited before
  // execution ever got here, so the check was removed from this copy
  // rather than kept as dead code.
  AN_CUDA_CHECK(cudaMalloc(&dtp, ngrid * niv * 2 * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&drp, ngrid * niv * 2 * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&ds, ce -> grid[1] * ce -> grid[2] * 5 * sizeof(double)));

  s = new double[5];
  a = new double[4];
  c = new double[2];
  x = new double[2];
  irho = new double[niv];
  for (int iv = 0; iv < niv; ++iv) {
    irho[iv] = 1.0 / sv -> rhov[iv];
  }
  initialized = true;
}


void AN2 :: calculate (double * dt, double * dtr) {
  __global__ void newdt0(double *, const double * __restrict__,
			 double *, double *);
  __global__ void newdt(double *, const double * __restrict__,
		       double *, double *,
		      double s1, double s2, double m, int niv, int biv);

  if (count > 0) {
    --count;
    for (int iv = 0; iv < niv; ++iv) {
      int biv = iv + binary * niv;
      newdt0 <<< g, b >>> (dt + (iv * ngrid), dtr + (iv * ngrid),
			    dtp + (biv * ngrid), drp + (biv * ngrid));
      AN_CUDA_CHECK(cudaGetLastError());
    }
  } else {
    cal_theta(dt, dtr);
    a[2] = a[1];
    leg(a, x, c, 2, 2);
    s1 = x[0] + mp * ((binary == 0)? 0 : 1);
    s2 = x[1] + mp * ((binary == 0)? 1 : 0);
    for (int iv = 0; iv < niv; ++iv) {
      newdt <<< g, b >>> (dt + (iv * ngrid), dtr + (iv * ngrid),
			   dtp + (iv * ngrid), drp + (iv * ngrid),
			   s1, s2, m, niv * ngrid, binary * niv * ngrid);
      AN_CUDA_CHECK(cudaGetLastError());
    }
  }
  binary = (binary == 0)? 1 : 0;
}


void AN2 :: cal_theta (double * dt, double * dtr) {
  __global__ void theta30(double *, const double * __restrict__,
                          const double * __restrict__, int);
  c[0] = c[1] = a[0] = a[1] = a[3] = 0.0;
  for (int iv = 0; iv < niv; ++iv) {
    theta30 <<< g, b, b.x * 5 * sizeof(double) >>>
      (ds, dtr + (iv * ngrid), drp + (iv * ngrid), niv * ngrid);
    AN_CUDA_CHECK(cudaGetLastError());
    thrust::device_ptr<double> ds_ptr(ds);
    for (int i = 0; i < 5; ++i) {
      s[i] = thrust::reduce(ds_ptr + i * g.x * g.y,
			    ds_ptr + (i + 1) * g.x * g.y);
    }
    c[0] += s[0] * irho[iv];
    c[1] += s[1] * irho[iv];
    a[0] += s[2] * irho[iv];
    a[1] += s[3] * irho[iv];
    a[3] += s[4] * irho[iv];
  }
}


static void leg(double * a, double * x, double * b, int n, int nn) {
  int i, j, k, k1;
  double s, p, q, r;
  double *ai, *ak;
  for (k = 0, ak = a; k < n -1; ++k, ak += nn) {
    k1 = k + 1;
    p = ak[k];
    for (j = k1; j < n; ++j)
      ak[j] /= p;
    r = b[k] /= p;
    for (i = k1, ai = ak + nn; i < n; ++i, ai += nn) {
      q = ai[k];
      for (j = k1; j < n; ++j)
        ai[j] -= q * ak[j];
      b[i] -= q * r;
    }
  }
  x[n - 1] = b[n - 1] / ak[n - 1];
  for (k = n - 2, ak = a + nn * (n - 2); k >= 0; --k, ak -= nn) {
    k1 = k + 1;
    s = b[k];
    for (j = k1; j < n; ++j)
      s -= ak[j] * x[j];
    x[k] = s;
  }
}
