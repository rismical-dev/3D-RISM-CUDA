#include <iostream>
#include <fstream>
#include <string>

#include "solvent.h"
#include "spline.h"
#include "alloc.h"
#include "cuda_check.h"

void Solvent :: spline2 (std::vector <double> & ga, int * & indga,
			int nga, int ngrid) {
  if (ga[nga - 1] > ttab2[ntab - 1]) {
    std::cout << "insufficient maximal T tabulated" << std::endl;
    exit (1);
  }

  std::vector <double *> chsa;
  alloc2D(chsa, natv, nga);
  alloc2D(wfka, natv, nga);

  int ntabb = 0;
  int ntabe = ntab2 - 1;
  int np = ntabe - ntabb + 1;

  double * x = new double[np];
  double * y = new double[np];
  std::vector <double *> coe;
  alloc2D(coe, 3, np);

  for (int n = 0; n < np; ++n) {
    x[n] = ttab2[ntabb + n] ;
  }

  for (int iv = 0; iv < natv; ++iv) {
    for (int n = 0; n < np; ++n) {
      y[n] = chs[iv][ntabb + n] ;
    }
    ::spline(x, y, np, coe) ;
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      chsa[iv][i] = splint(x, y, coe, np, ga[i]);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
    for (int n = 0; n < np; ++n) {
      y[n] = wfk[iv][ntabb + n] ;
    }
    ::spline(x, y, np, coe) ;
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      wfka[iv][i] = splint(x, y, coe, np, ga[i]);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      cvva[iv][iv][i] -= chsa[iv][i] * rhov[0] / rhov[iv];
    }
  }

  AN_CUDA_CHECK(cudaMalloc(&dc, nga * natv * natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&dw, nga * natv * sizeof(double)));

  // Synchronous, since cvva is deallocated (dealloc3D) right below: an
  // async copy here would let that free race the H2D transfer (same class
  // of bug fixed for siguv/epsuv in cal_LJ.cu).
  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
      AN_CUDA_CHECK(cudaMemcpy(dc + (iv1 * nga) + (iv2 * natv * nga),
      		      cvva[iv2][iv1], nga * sizeof(double),
      		      cudaMemcpyDefault));
    }
  }

  // wfka is a member (kept alive after this function returns, unlike
  // cvva/chsa/coe above), so the async copy here is safe: nothing frees
  // wfka before whatever later kernel actually consumes dw.
  for (int iv = 0; iv < natv; ++iv) {
    AN_CUDA_CHECK(cudaMemcpyAsync(dw + (iv * nga), wfka[iv], nga * sizeof(double),
		    cudaMemcpyDefault));
  }

  // rhov is also a member kept alive after this function returns, so this
  // async copy is likewise safe.
  AN_CUDA_CHECK(cudaMalloc(&drho, natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMemcpyAsync(drho, rhov, natv * sizeof(double), cudaMemcpyDefault));

  dealloc3D(cvva);
  dealloc2D(coe);
  dealloc2D(chsa);
  delete[] x;
  delete[] y;
  delete[] ttab2;
}
