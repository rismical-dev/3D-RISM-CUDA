#include <thrust/device_vector.h>
#include "rism3d.h"

void RISM3D :: cal_qv (double * & hqv) {
  __global__ void qv(double *, const double2 * __restrict__,
  	     double);
  double * dqv;
  cudaMalloc(&dqv, ce -> ngrid * sizeof(double));
  cudaMemset(dqv, 0.0, ce -> ngrid * sizeof(double));

  for (int iv = 0; iv < sv -> natv; ++iv) {
    qv <<< g, b >>> (dqv, dguv + (iv * ce -> ngrid),
       	      	     sv -> rhov[iv] * sv -> qv[iv]);
  }

  cudaMemcpy(hqv, dqv, ce -> ngrid * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(dqv);
} 

__global__ void qv(double * dqv, const double2 * __restrict__ dguv,
	           double rq) {

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;

  dqv[ip] += dguv[ip].x * rq;
}
