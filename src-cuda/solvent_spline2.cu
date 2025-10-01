#include <iostream>
#include <fstream>
#include <string>

#include "solvent.h"

void Solvent :: spline2 (vector <double> & ga, int * & indga,
			int nga, int ngrid) {
  void alloc2D (vector <double *> &, int, int);
  void dealloc2D (vector <double *> &);
  void dealloc3D (vector < vector <double *> > &);

  void spline (double * &, double * &, int, vector <double *> &);
  double splint (double * &, double * &, vector <double *> &, int, double);

  if (ga[nga - 1] > ttab2[ntab - 1]) {
    cout << "insufficient maximal T tabulated" << endl;
    exit (1);
  }

  vector <double *> chsa;
  alloc2D(chsa, natv, nga);
  alloc2D(wfka, natv, nga);

  int ntabb = 0;
  for (int n = 0; n < ntab2; ++n) {
    if (ttab2[n] < ga[0]) 
      ntabb = n;
  }

  int ntabe = ntab2 - 1;
  for (int n = ntab2 - 1; n > 0; --n) {
    if (ttab2[n] > ga[nga - 1]) 
      ntabe = n;
  }

  int np = ntabe - ntabb + 1;

  double * x = new double[np];
  double * y = new double[np];
  vector <double *> coe;
  alloc2D(coe, 3, np);

  for (int n = 0; n < np; ++n) {
    x[n] = ttab2[ntabb + n] ;
  }

  for (int iv = 0; iv < natv; ++iv) {
    for (int n = 0; n < np; ++n) {
      y[n] = chs[iv][ntabb + n] ;
    }
    spline(x, y, np, coe) ;
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      chsa[iv][i] = splint(x, y, coe, np, ga[i]);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
    for (int n = 0; n < np; ++n) {
      y[n] = wfk[iv][ntabb + n] ;
    }
    spline(x, y, np, coe) ;
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

  cudaMalloc(&dc, nga * natv * natv * sizeof(double));
  cudaMalloc(&dw, nga * natv * sizeof(double));
  
  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
      cudaMemcpyAsync(dc + (iv1 * nga) + (iv2 * natv * nga),
      		      cvva[iv2][iv1], nga * sizeof(double),
      		      cudaMemcpyDefault);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
    cudaMemcpyAsync(dw + (iv * nga), wfka[iv], nga * sizeof(double),
		    cudaMemcpyDefault);
  }

  cudaMalloc(&drho, natv * sizeof(double));
  cudaMemcpyAsync(drho, rhov, natv * sizeof(double), cudaMemcpyDefault);

  dealloc3D(cvva);
  dealloc2D(coe);
  dealloc2D(chsa);
  delete[] x;
  delete[] y;
  delete[] ttab2;
}
