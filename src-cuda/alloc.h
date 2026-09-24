#ifndef ALLOC_H
#define ALLOC_H
#include <vector>
#include <complex>

// Manual "2D/3D array" helpers built on vector<double*> /
// vector<vector<double*>> / vector<complex<double>*>, used throughout the
// codebase (solvent_read.cc, solvent_spline.cu, solvent_spline2.cu,
// solvent.cu, iterate.cu, ...). Declared here instead of each caller
// locally forward-declaring them, so they can't drift out of sync with
// alloc2D.cc/alloc3D.cc/calloc2D.cc/dealloc2D.cc/dealloc3D.cc's actual
// signatures.
//
// alloc2D/alloc3D leave their elements uninitialized (plain new double[]),
// matching the original code. calloc2D's elements come out zero-valued
// only because std::complex<double>'s default constructor does that, not
// because of anything calloc2D itself does -- despite the name, it is not
// a zeroing allocation the way C's calloc() is.
//
// dealloc2D is overloaded for both row types calloc2D/alloc2D can produce,
// since a vector<complex<double>*> can't bind to the vector<double*>
// overload.
void alloc2D (std::vector <double *> & array2, int r1, int r2);
void alloc3D (std::vector <std::vector <double *> > & array3,
              int r1, int r2, int r3);
void calloc2D (std::vector <std::complex <double> *> & array2,
               int r1, int r2);
void dealloc2D (std::vector <double *> & array2);
void dealloc2D (std::vector <std::complex <double> *> & array2);
void dealloc3D (std::vector <std::vector <double *> > & array3);

#endif  // ALLOC_H
