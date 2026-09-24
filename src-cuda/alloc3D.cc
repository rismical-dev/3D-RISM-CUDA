#include "alloc.h"
using namespace std;

void alloc3D (vector < vector <double* > > & array3, int r1, int r2, int r3) {
  for (int r = 0 ; r < r1 ; ++r) {
    vector <double *> temp2;
    // Inner loop variable renamed from the original's "r" (which shadowed
    // the outer loop's "r") to "c" so the two can't be confused.
    for (int c = 0 ; c < r2 ; ++c) {
      double * temp1 = new double[r3];
      temp2 . push_back (temp1) ;
    }
    array3 . push_back (temp2) ;
  }
}
