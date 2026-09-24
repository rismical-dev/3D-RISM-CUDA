#include <thrust/device_vector.h>
#include "rism3d.h"
#include "grid_constants.h"

// See cal_ad1.cu's comment: renamed from "du" (which shadowed the
// unrelated RISM3D device member "du", the LJ potential grid) to "ad"
// (atomic decomposition), matching the array's name at its call site in
// output.cu.
void RISM3D :: cal_ad2(double * & ad) {
  __global__ void ad2(double * ds, double2 * dguv, double * dsig,
	              double * deps,  double3 * dr, double * qu,
		      double qv, int natu, int iv, int iu, double lambda);

  int ng = ce -> ngrid;

  set_grid_constants(ce);

  double * ds;
  AN_CUDA_CHECK(cudaMalloc(&ds, g.x * g.y * sizeof(double)));

#pragma omp parallel for
  for (int iu = 0; iu < su -> num; ++iu) {
    ad[iu] = 0.0;
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int iu = 0; iu < su -> num; ++iu) {
      ad2 <<< g, b, b.x * sizeof(double) >>>
	(ds, dguv + (iv * ng), dsig, deps, su -> dr, su -> dq, sv -> qv[iv],
	 su -> num, iv, iu, lambda);
      AN_CUDA_CHECK(cudaGetLastError());

      thrust::device_ptr<double> ds_ptr(ds);
      double s = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
      ad[iu] += s * sv -> rhov[iv];
    }
  }
  AN_CUDA_CHECK(cudaFree(ds));
}

__global__ void ad2(double * ds, double2 * dguv, double * dsig,
		    double * deps,  double3 * dr, double * qu,
                    double qv, int natu, int iv, int iu, double lambda) {
  extern __shared__ double sdata[];
  const double cc = hartree * bohr * avogadoro;

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  int iuv = iu + iv * natu;

  double dx = ((int)threadIdx.x - grid.x / 2) * dv.x - dr[iu].x;
  double dy = ((int)blockIdx.x - grid.y / 2) * dv.y - dr[iu].y;
  double dz = ((int)blockIdx.y - grid.z / 2) * dv.z - dr[iu].z;
  double r2 = dx * dx + dy * dy + dz * dz;
  double r1 = sqrt(r2);

  if (r1 < dsig[iuv] * 0.5) {
    sdata[threadIdx.x] = 0.0;
  } else {
    double uco = qu[iu] * qv / r1 * cc * dguv[ip].x;
    sdata[threadIdx.x] = uco;
  }
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
  }
}
