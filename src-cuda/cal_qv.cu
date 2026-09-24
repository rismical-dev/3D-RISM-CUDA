#include <thrust/device_vector.h>
#include "rism3d.h"
#include "cuda_check.h"

// Same indexing formula as fft3d.cu / anderson_cuda.cu / initialize_g.cu's
// gridIndex() helper. Kept as its own static (internal-linkage) copy here
// rather than sharing a cross-TU declaration, for the same reason those
// files each have their own copy instead of one shared one: it's a
// one-line __device__ helper used only by this file's own kernel.
static __device__ __forceinline__ unsigned int gridIndex() {
  return threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
}

void RISM3D :: cal_qv (double * & hqv) {
  __global__ void qv(double *, const double2 * __restrict__,
  	     double);
  double * dqv;
  AN_CUDA_CHECK(cudaMalloc(&dqv, ce -> ngrid * sizeof(double)));
  AN_CUDA_CHECK(cudaMemset(dqv, 0, ce -> ngrid * sizeof(double)));

  for (int iv = 0; iv < sv -> natv; ++iv) {
    qv <<< g, b >>> (dqv, dguv + (iv * ce -> ngrid),
       	      	     sv -> rhov[iv] * sv -> qv[iv]);
    AN_CUDA_CHECK(cudaGetLastError());
  }

  AN_CUDA_CHECK(cudaMemcpy(hqv, dqv, ce -> ngrid * sizeof(double), cudaMemcpyDeviceToHost));
  AN_CUDA_CHECK(cudaFree(dqv));
}

__global__ void qv(double * dqv, const double2 * __restrict__ dguv,
	           double rq) {

  unsigned int ip = gridIndex();

  dqv[ip] += dguv[ip].x * rq;
}
