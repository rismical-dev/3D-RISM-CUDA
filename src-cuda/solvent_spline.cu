#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

#include "solvent.h"

void Solvent :: spline (vector <double> & ga, int * & indga,
			int nga, int ngrid, bool rmdft) {
  void alloc2D (vector <double *> &, int, int);
  void alloc3D (vector <vector <double *> > &, int, int, int);
  void dealloc2D (vector <double *> &);
  void dealloc3D (vector < vector <double *> > &);

  void spline (double * &, double * &, int, vector <double *> &);
  double splint (double * &, double * &, vector <double *> &, int, double);

  if (ga[nga - 1] > ttab[ntab - 1]) {
    cout << "insufficient maximal T tabulated" << endl;
    exit (1);
  }

  vector <vector <double *> > xvva2;
  alloc3D (xvva2, natv, natv, nga);
  if (rmdft) alloc3D (cvva, natv, natv, nga);

  int ntabb = 0;
  int ntabe = ntab - 1;
  int np = ntabe - ntabb + 1;

  double * x = new double[np];
  double * y = new double[np];
  vector <double *> coe;
  alloc2D(coe, 3, np);

  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
      for (int n = 0; n < np; ++n) {
	x[n] = ttab[ntabb + n] ;
	y[n] = xvv[iv2][iv1][ntabb + n] ;
      }
      spline(x, y, np, coe) ;
#pragma omp parallel for
      for (int i = 0; i < nga; ++i) {
	xvva2[iv2][iv1][i] = splint(x, y, coe, np, ga[i]);
      }
    }
  }

  if (rmdft) {
    for (int iv2 = 0; iv2 < natv; ++iv2) {
      for (int iv1 = 0; iv1 < natv; ++iv1) {
        for (int n = 0; n < np; ++n) {
           y[n] = cvv[iv2][iv1][ntabb + n] ;
        }
        spline(x, y, np, coe) ;
#pragma omp parallel for
        for (int i = 0; i < nga; ++i) {
          cvva[iv2][iv1][i] = splint(x, y, coe, np, ga[i]);
	}
      }
    }
  }

  alloc3D(xvva, natv, natv, ngrid);
  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
#pragma omp parallel for
      for (int ig = 0; ig < ngrid; ++ig) {
	int iga = indga[ig];
	xvva[iv2][iv1][ig] = xvva2[iv2][iv1][iga];
      }
    }
  }

  cudaMalloc(&dx, ngrid * natv * natv * sizeof(double));
  
  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
      cudaMemcpyAsync(dx + (iv1 * ngrid) + (iv2 * natv * ngrid),
      		      xvva[iv2][iv1], ngrid * sizeof(double),
      		      cudaMemcpyDefault);
    }
  }

  dealloc3D(xvva2);
  dealloc3D(xvv);
  dealloc2D(coe);
  delete[] x, y;
  delete[] ttab;
}
