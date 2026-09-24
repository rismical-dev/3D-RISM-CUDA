#include "alloc.h"
using namespace std;

void dealloc2D (vector <double *> & array2) {
  for (vector<double *>::size_type r = 0; r < array2.size(); ++r) {
    delete[] array2[r];
  }
  array2.clear();
}

// Overload for the complex<double>* rows calloc2D() produces (e.g.
// RISM3D's guv/huv in iterate.cu). The vector<double*> overload above
// can't accept these since the element type doesn't match, which is why
// this overload didn't exist before -- and why calloc2D-allocated arrays
// had no matching deallocation function anywhere in the codebase.
void dealloc2D (vector <complex <double> *> & array2) {
  for (vector<complex <double> *>::size_type r = 0; r < array2.size(); ++r) {
    delete[] array2[r];
  }
  array2.clear();
}
