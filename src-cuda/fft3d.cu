#include <iostream>
#include "fft3d.h"
#include "cuda_check.h"
#include "cufft_check.h"

// Every kernel in this file addresses the same flattened 3D grid the same
// way: b.x/g.x/g.y (set up in FFT3D::initialize) give threadIdx.x the
// fastest-varying axis and blockIdx.x/blockIdx.y the other two. Factored
// out so the indexing scheme only has to be read (and fixed, if it ever
// needs to change) in one place instead of six.
__device__ __forceinline__ unsigned int gridIndex () {
  return threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
}

FFT3D :: ~FFT3D () {
  if (!initialized) return;
  // Best-effort cleanup: a destructor shouldn't abort the whole run over a
  // teardown failure (unlike AN_CUDA_CHECK/AN_CUFFT_CHECK in initialize()/
  // execute(), where a failure means the FFT itself can no longer be
  // trusted), so this reports problems instead of exit()-ing on them.
  cufftResult fr = cufftDestroy(plan);
  if (fr != CUFFT_SUCCESS) {
    fprintf(stderr, "cuFFT warning at %s:%d: cufftDestroy failed: %s (%d)\n",
            __FILE__, __LINE__, anCufftErrorString(fr), (int)fr);
  }
  cudaError_t er = cudaFree(dkf);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dkf) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
  er = cudaFree(dir);
  if (er != cudaSuccess) {
    fprintf(stderr, "CUDA/HIP warning at %s:%d: cudaFree(dir) failed: %s\n",
            __FILE__, __LINE__, cudaGetErrorString(er));
  }
}

void FFT3D :: initialize (Cell * ce) {
  __global__ void set_kf(double2 * dkf, int nx, int ny, int nz);
  __global__ void set_ir(int * dir);

  ngrid = ce -> ngrid;
  volf = ce -> dv;
  volb = 1.0 / ce -> volume;
  b.x = ce -> grid[0];
  g.x = ce -> grid[1];
  g.y = ce -> grid[2];
  AN_CUDA_CHECK(cudaMalloc(&dkf, ngrid * sizeof(double2)));
  AN_CUDA_CHECK(cudaMalloc(&dir, ngrid * sizeof(int)));
  AN_CUFFT_CHECK(cufftPlan3d(&plan, ce -> grid[0], ce -> grid[1], ce -> grid[2], CUFFT_Z2Z));
  set_kf <<< g, b >>> (dkf, ce -> grid[0], ce -> grid[1], ce -> grid[2]);
  AN_CUDA_CHECK(cudaGetLastError());
  set_ir <<< g, b >>> (dir);
  AN_CUDA_CHECK(cudaGetLastError());
  initialized = true;
}

void FFT3D :: execute (double2 * da, int key) {
  __global__ void timeirvol(double2 *, const int * __restrict__, double);
  __global__ void timekf(double2 *, const double2 * __restrict__,
  			 const int * __restrict__);
  __global__ void timekb(double2 *, const double2 * __restrict__,
  			 const int * __restrict__, double);

  if (key == FORWARD) {
    timekf <<< g, b >>> (da, dkf, dir);
    AN_CUDA_CHECK(cudaGetLastError());
    AN_CUFFT_CHECK(cufftExecZ2Z(plan, da, da, CUFFT_FORWARD));
    timeirvol <<< g, b >>> (da, dir, volf);
    AN_CUDA_CHECK(cudaGetLastError());
  } else {
    // Anything other than FORWARD is treated as BACKWARD, exactly as the
    // original "if (key == -1) {...} else {...}" did.
    timeirvol <<< g, b >>> (da, dir, 1.0);
    AN_CUDA_CHECK(cudaGetLastError());
    AN_CUFFT_CHECK(cufftExecZ2Z(plan, da, da, CUFFT_INVERSE));
    timekb <<< g, b >>> (da, dkf, dir, volb);
    AN_CUDA_CHECK(cudaGetLastError());
  }
}

__global__ void timekf(double2 * da, const double2 * __restrict__ dkf,
		       const int * __restrict__ dir) {
  unsigned int ip = gridIndex();
  double tmpr = da[ip].x * dir[ip];
  double tmpi = da[ip].y * dir[ip];
  da[ip].x = tmpr * dkf[ip].x - tmpi * dkf[ip].y;
  da[ip].y = tmpi * dkf[ip].x + tmpr * dkf[ip].y;
}

__global__ void timekb(double2 * da, const double2 * __restrict__ dkf,
		       const int * __restrict__ dir, double vol) {
  unsigned int ip = gridIndex();
  double tmpr = da[ip].x * dir[ip] * vol;
  double tmpi = da[ip].y * dir[ip] * vol;
  da[ip].x = tmpr * dkf[ip].x + tmpi * dkf[ip].y;
  da[ip].y = tmpi * dkf[ip].x - tmpr * dkf[ip].y;
}

__global__ void timeirvol(double2 * da, const int * __restrict__ dir, double vol) {
  unsigned int ip = gridIndex();
  da[ip].x *= dir[ip] * vol;
  da[ip].y *= dir[ip] * vol;
}

__global__ void set_kf(double2 * dkf, int nx, int ny, int nz) {
  unsigned int ip = gridIndex();
  int x = threadIdx.x - nx / 2;
  int y = blockIdx.x - ny / 2;
  int z = blockIdx.y - nz / 2;
  double dkx = M_PI / nx;
  double dky = M_PI / ny;
  double dkz = M_PI / nz;
  double dkr = dkx * x + dky * y + dkz * z;
  dkf[ip].x = cos(dkr);
  dkf[ip].y = - sin(dkr);
}

__global__ void set_ir(int * dir) {
  unsigned int ip = gridIndex();
  if ((threadIdx.x + blockIdx.x + blockIdx.y) % 2 == 0) {
    dir[ip] = 1;
  } else {
    dir[ip] = -1;
  }
}
