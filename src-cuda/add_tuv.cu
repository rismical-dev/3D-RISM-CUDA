#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "cuda_check.h"

void RISM3D :: add_tuv (double cuf) {
  __global__ void add_dt(double * dt, const double * __restrict__ fr, double q);

  for (int iv = 0; iv < sv -> natv; ++iv) {
    add_dt <<< g, b >>> (dt + (iv * ce -> ngrid), dfr, sv -> qv[iv] * cuf);
    AN_CUDA_CHECK(cudaGetLastError());
  }
}

__global__ void add_dt(double * dt, const double * __restrict__ fr, double q) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  dt[ip] += q * fr[ip];
}
