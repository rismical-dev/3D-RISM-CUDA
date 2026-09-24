#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include "rism3d.h"

template <typename T> struct square {
  __host__ __device__ T operator()(const T &x) const {
    return x * x;
  }
};


double RISM3D :: cal_rms () {
  square<double> uop;
  thrust::plus<double> bop;
  thrust::device_ptr<double> dtr_ptr(dtr);

  // size_t, not int * int: ngrid * natv can exceed INT_MAX for large
  // enough grids, which would otherwise overflow both the iterator bound
  // below and the divisor in the rms normalization.
  size_t n = static_cast<size_t>(ce -> ngrid) * sv -> natv;

  double rms = thrust::transform_reduce(dtr_ptr, dtr_ptr
				+ n, uop, 0.0, bop);
  rms = sqrt (rms / n);
  return rms;
}
