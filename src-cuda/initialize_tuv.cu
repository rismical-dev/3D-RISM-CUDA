#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "cuda_check.h"

// Same indexing formula as cal_qv.cu / cal_grad.cu / cal_euv.cu /
// fft3d.cu / anderson_cuda.cu / initialize_g.cu's gridIndex() helper. Kept
// as its own static (internal-linkage) copy here rather than sharing a
// cross-TU declaration, for the same reason those files each have their
// own copy.
static __device__ __forceinline__ unsigned int gridIndex() {
  return threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
}

void RISM3D :: initialize_tuv () {
  __global__ void init_tuv(double * dt, double * fr, double q);

  std::cout << "synthesizing initial estimate for Tuv ..." << std::endl ;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    init_tuv <<< g, b >>> (dt + (iv * ce -> ngrid), dfr, sv -> qv[iv]);
    AN_CUDA_CHECK(cudaGetLastError());
  }
}

__global__ void init_tuv(double * dt, double * fr, double q) {
  unsigned int ip = gridIndex();
  dt[ip] = q * fr[ip];
}
