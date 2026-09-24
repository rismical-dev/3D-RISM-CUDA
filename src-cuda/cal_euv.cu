#include <thrust/device_vector.h>
#include "rism3d.h"
#include "cuda_check.h"
#include "grid_constants.h"

// Same indexing formula as cal_qv.cu / cal_grad.cu / fft3d.cu /
// anderson_cuda.cu / initialize_g.cu's gridIndex() helper. Kept as its own
// static (internal-linkage) copy here rather than sharing a cross-TU
// declaration, for the same reason those files each have their own copy.
static __device__ __forceinline__ unsigned int gridIndex() {
  return threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
}

void RISM3D :: cal_euv (double * & e) {
  __global__ void euv(double * ds, double2 * dhuv, double * dsig,
                      double * deps,  double3 * dr, double * qu, double gv,
		      int natu, int iv, int iu);

  int ng = ce -> ngrid;

  set_grid_constants(ce);
  double * ds2;
  AN_CUDA_CHECK(cudaMalloc(&ds2, g.x * g.y * 2 * sizeof(double)));

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int iu = 0; iu < su -> num; ++iu) {
      euv <<< g, b, b.x * 2 * sizeof(double) >>>
        (ds2, dguv + iv * ng, dsig, deps, su -> dr, su -> dq, sv -> qv[iv],
         su -> num, iv, iu);
      AN_CUDA_CHECK(cudaGetLastError());
      thrust::device_ptr<double> ds2_ptr(ds2);
      for (int i = 0; i < 2; ++i) {
        double s = thrust::reduce(ds2_ptr + (g.x * g.y) * i,
                                ds2_ptr + (g.x * g.y) * (i + 1));
        e[(sv -> natv * su -> num) * i + su -> num * iv + iu] =
  	  s * sv -> rhov[iv];
      }
    }
  }
  AN_CUDA_CHECK(cudaFree(ds2));
}

__global__ void euv(double * ds, double2 * dguv, double * dsig,
                    double * deps,  double3 * dr, double * qu,
		    double qv, int natu, int iv, int iu) {
  extern __shared__ double sdata[];
  const double cc = hartree * bohr * avogadoro;

  unsigned int ip = gridIndex();
  int iuv = iu + iv * natu;

  double dx = ((int)threadIdx.x - grid.x / 2) * dv.x - dr[iu].x;
  double dy = ((int)blockIdx.x - grid.y / 2) * dv.y - dr[iu].y;
  double dz = ((int)blockIdx.y - grid.z / 2) * dv.z - dr[iu].z;
  double r2 = dx * dx + dy * dy + dz * dz;
  double r1 = sqrt(r2);

  if (r1 < dsig[iuv] * 0.5) {
    sdata[threadIdx.x] = 0.0;
    sdata[threadIdx.x + blockDim.x] = 0.0;
  } else {
    double rs2i = dsig[iuv] * dsig[iuv] / r2;
    double rs6i = rs2i * rs2i * rs2i;
    double ulj = deps[iuv] * 4.0 * rs6i * ( rs6i - 1.0) * dguv[ip].x;
    double uco = qu[iu] * qv / r1 * cc * dguv[ip].x;
    sdata[threadIdx.x] = ulj;
    sdata[threadIdx.x + blockDim.x] = uco;
  }
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x + blockDim.x] += sdata[threadIdx.x + blockDim.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
    ds[blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y] =
      sdata[blockDim.x];
  }
}
