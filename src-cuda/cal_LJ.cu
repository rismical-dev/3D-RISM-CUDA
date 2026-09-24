#include <iostream>
#include "rism3d.h"
#include "cuda_check.h"

void RISM3D :: cal_LJ() {
  __global__ void LJ(double * du, const double * __restrict__ dsig,
		     const double * __restrict__ deps,
		     const double3 * __restrict__ dru,
		     double cut2, double ikbt, double bx, double by, double bz,
		     int nx, int ny, int nz, int natu, int iv);

  const double cut = 2.0e-3;
  const double cut2 = cut * cut;

  std::cout << "tabulating solute Lennard-Jones potential ..." << std::endl;

  AN_CUDA_CHECK(cudaMalloc(&du, ce -> ngrid * sv -> natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&dsig, su -> num * sv -> natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&deps, su -> num * sv -> natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMemset(du, 0, ce -> ngrid * sv -> natv * sizeof(double)));

  // Local staging buffers: only ever used here to build the host-side
  // sig/eps combination rules before copying them to the device, and
  // nowhere else in the codebase (confirmed by grep). Freed at the end of
  // this function instead of being kept around as unused RISM3D members
  // for the rest of the run.
  double * siguv = new double[su -> num * sv -> natv];
  double * epsuv = new double[su -> num * sv -> natv];

  double lambda1 = 1.0;
  if (adswitch == 1) lambda1 = lambda;

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int iu = 0; iu < su -> num; ++iu) {
      int ip = iu + su -> num * iv;
      siguv[ip] = (su -> sig[iu] + sv -> sigv[iv]) * 0.5 * lambda1;
      epsuv[ip] = sqrt (su -> eps[iu] * sv -> epsv[iv]);
    }
  }

  // Synchronous, unlike the fire-and-forget H2D copies elsewhere in this
  // codebase that are immediately followed only by device-side kernel
  // launches (safe under stream ordering without a host-side wait): here
  // siguv/epsuv are deleted right below, so the copy must actually have
  // finished reading them before that delete runs.
  AN_CUDA_CHECK(cudaMemcpy(dsig, siguv, su -> num * sv -> natv * sizeof(double),
			  cudaMemcpyDefault));
  AN_CUDA_CHECK(cudaMemcpy(deps, epsuv, su -> num * sv -> natv * sizeof(double),
			  cudaMemcpyDefault));
  delete[] siguv;
  delete[] epsuv;

  double iKbT = 1.0 / (avogadoro * boltzmann * sv -> temper);
  for (int iv = 0; iv < sv -> natv; ++iv) {
    LJ <<< g, b >>> (du + (iv * ce -> ngrid), dsig, deps, su -> dr,
		      cut2, iKbT, ce -> dr[0], ce -> dr[1], ce -> dr[2],
		      ce -> grid[0], ce -> grid[1], ce -> grid[2],
		      su -> num, iv);
    AN_CUDA_CHECK(cudaGetLastError());
  }
}

__global__ void LJ(double * du, const double * __restrict__ dsig,
		   const double * __restrict__ deps,
		   const double3 * __restrict__ dru,
                   double cut2, double ikbt, double bx, double by, double bz,
                   int nx, int ny, int nz, int natu, int iv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rx = ((int)threadIdx.x - nx / 2) * bx;
  double ry = ((int)blockIdx.x - ny / 2) * by;
  double rz = ((int)blockIdx.y - nz / 2) * bz;
  for (int iu = 0; iu < natu; ++iu) {
    int iuv = iu + natu * iv;
    double dx = rx - dru[iu].x;
    double dy = ry - dru[iu].y;
    double dz = rz - dru[iu].z;
    double r2 = dx * dx + dy * dy + dz * dz ;

    double rs2 = r2 / (dsig[iuv] * dsig[iuv]);
    if (rs2 < cut2) rs2 = cut2;
    double irs6 = 1.0 / (rs2 * rs2 * rs2);
    du[ip] += deps[iuv] * 4.0 * irs6 * (irs6 - 1.0);
  }
  du[ip] *= ikbt;
}
