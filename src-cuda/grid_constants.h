#ifndef GRID_CONSTANTS_H
#define GRID_CONSTANTS_H

#include "cell.h"
#include "cuda_check.h"

// Shared by the per-atom decomposition kernels (cal_ad1.cu, cal_ad2.cu,
// cal_grad.cu, cal_euv.cu): the cell's grid spacing and dimensions, needed
// in constant memory by kernels that compute a grid point's position from
// thread/block indices.
//
// Marked static (rather than plain __constant__) to force internal linkage
// explicitly: this codebase doesn't build with -rdc=true (relocatable
// device code), so each .cu file that includes this header already gets
// its own private copy with no cross-file linker conflict, but static
// keeps that true even if -rdc were ever turned on later.
static __constant__ double3 dv;
static __constant__ int3 grid;

// static (not just inline): each .cu file that includes this header must
// get its OWN private copy of this function bound to its OWN private
// static dv/grid symbols above. A plain "inline" here would still be
// legal C++, but a toolchain that treats "inline" as license to fold
// identical-looking definitions from different translation units into a
// single kept copy (COMDAT-style deduplication) can end up keeping only
// one file's copy and routing every other file's call to it -- which
// then writes into that one file's dv/grid instead of the caller's own,
// leaving the caller's kernels reading never-initialized (zero) grid
// spacing/dimensions. That was diagnosed as the cause of a real
// regression in cal_euv.cu's SE output once a 4th file (after
// cal_ad1.cu/cal_ad2.cu/cal_grad.cu) started sharing this header, and
// "static" is the fix: it forces true internal linkage, so the compiler
// can never merge this definition across translation units.
static inline void set_grid_constants(Cell * ce) {
  AN_CUDA_CHECK(cudaMemcpyToSymbol(dv, ce -> dr, sizeof(double3)));
  AN_CUDA_CHECK(cudaMemcpyToSymbol(grid, ce -> grid, sizeof(int3)));
}

#endif  // GRID_CONSTANTS_H
