#ifndef Cell_H
#define Cell_H

class Cell {
public:
  Cell() {
    box = new double[3];
    shift = new double[3];
    dr = new double[3];
    grid = new int[3];
  }
  ~Cell() {delete[] box; delete[] shift; delete[] dr; delete[] grid;}

  // Cell owns four heap-allocated arrays. Copying it would give two
  // objects the same four pointers, and whichever is destroyed first
  // would free memory the other still thinks it owns. Nothing in the
  // codebase copies a Cell today (RISM3D only ever holds one through a
  // pointer), so disabling copy costs nothing and rules that bug class
  // out entirely.
  Cell (const Cell &) = delete;
  Cell & operator= (const Cell &) = delete;

  void setup();
  double * box;
  double * shift;
  double * dr;
  int * grid;
  double volume, dv;
  int ngrid;
};

#endif // Cell_H
