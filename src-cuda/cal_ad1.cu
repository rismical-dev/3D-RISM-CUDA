#include <thrust/device_vector.h>
#include "rism3d.h"
#include "grid_constants.h"

// The output parameter is named "ad" (atomic decomposition), matching the
// array's name at its call site in output.cu, instead of "du" -- "du" also
// happens to be the name of an unrelated RISM3D device member (the LJ
// potential grid built in cal_LJ.cu). The old name was never actually
// wrong (a function parameter always shadows a same-named member, so this
// function did read/write its own local array correctly either way), but
// it was easy to misread as touching that member.
void RISM3D :: cal_ad1(double * & ad) {
  __global__ void ad1(double * ds, double2 * dguv, double * dsig,
		       double * deps,  double3 * dr,
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
      ad1 <<< g, b, b.x * sizeof(double) >>>
	(ds, dguv + (iv * ng), dsig, deps, su -> dr, sv -> qv[iv],
	 su -> num, iv, iu, lambda);
      AN_CUDA_CHECK(cudaGetLastError());

      thrust::device_ptr<double> ds_ptr(ds);
      double s = thrust::reduce(ds_ptr, ds_ptr + g.x * g.y);
      ad[iu] += s * sv -> rhov[iv];
    }
  }
  AN_CUDA_CHECK(cudaFree(ds));
}

// qu (the solute charges) was a parameter here but unused in the kernel
// body -- this is the LJ-only decomposition, so it needs no charge data.
// It's kept out of the signature now instead of being carried along
// unused just to mirror cal_ad2's ad2() call, whose Coulomb decomposition
// does need it.
__global__ void ad1(double * ds, double2 * dguv, double * dsig,
		   double * deps,  double3 * dr,
                   double qv, int natu, int iv, int iu,
 		   double lambda) {
  extern __shared__ double sdata[];

  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  int iuv = iu + iv * natu;

  double dx = ((int)threadIdx.x - grid.x / 2) * dv.x - dr[iu].x;
  double dy = ((int)blockIdx.x - grid.y / 2) * dv.y - dr[iu].y;
  double dz = ((int)blockIdx.y - grid.z / 2) * dv.z - dr[iu].z;
  double r2 = dx * dx + dy * dy + dz * dz;
  double r1 = sqrt(r2);

  if (r1 < dsig[iuv] * lambda * 0.5) {
    sdata[threadIdx.x] = 0.0;
  } else {
    double rs2i = dsig[iuv] * dsig[iuv] / r2;
    double rs6i = rs2i * rs2i * rs2i;
    double ulj = deps[iuv] * 24.0 * rs6i * (2.0 * rs6i - 1.0) / lambda
      * dguv[ip].x;
    sdata[threadIdx.x] = ulj;
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
