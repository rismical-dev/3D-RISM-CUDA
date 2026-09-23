// Every kernel in this file addresses the same flattened 3D grid the same
// way (see AN2::initialize: b.x/g.x/g.y give threadIdx.x the fastest-
// varying axis and blockIdx.x/blockIdx.y the other two). Factored out so
// the indexing scheme only has to be read (and fixed, if it ever needs to
// change) in one place. static (internal linkage) because fft3d.cu
// defines its own same-named helper for the same reason -- keeping both
// static avoids a duplicate-symbol clash if these ever get compiled with
// device-code linking (-rdc=true) enabled.
static __device__ __forceinline__ unsigned int gridIndex () {
  return threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
}

__global__ void newdt0 (double * dt, const double * __restrict__ dtr,
                        double * dtp, double * drp) {
    unsigned int ip = gridIndex();
    dtp[ip] = dt[ip];
    drp[ip] = dtr[ip];
    dt[ip] += dtr[ip];
}


__global__ void newdt(double * dt, const double * __restrict__ dtr,
		      double * dtp, double * drp,
		      double s1, double s2, double m, int niv, int biv) {
  unsigned int ip = gridIndex();
  double u = dt[ip] + s1 * (dtp[ip] - dt[ip]) + s2 * (dtp[ip + niv] - dt[ip]);
  double v = dt[ip] + dtr[ip]
    + s1 * (dtp[ip] + drp[ip] - dt[ip] - dtr[ip])
    + s2 * (dtp[ip + niv] + drp[ip + niv] - dt[ip] - dtr[ip]);
  dtp[ip + biv] = dt[ip];
  drp[ip + biv] = dtr[ip];
  dt[ip] = u + m * (v - u);
}


__global__ void theta30(double * ds, const double * __restrict__ dtr,
			const double * __restrict__ drp, int niv) {
  extern __shared__ double sdata[];

  unsigned int ip = gridIndex();

  double t1 = dtr[ip] - drp[ip];
  double t2 = dtr[ip] - drp[ip + niv];
  sdata[threadIdx.x] = dtr[ip] * t1;
  sdata[threadIdx.x + blockDim.x] = dtr[ip] * t2;
  sdata[threadIdx.x + blockDim.x * 2] = t1 * t1;
  sdata[threadIdx.x + blockDim.x * 3] = t1 * t2;
  sdata[threadIdx.x + blockDim.x * 4] = t2 * t2;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      sdata[threadIdx.x] += sdata[threadIdx.x + s];
      sdata[threadIdx.x + blockDim.x] += sdata[threadIdx.x + blockDim.x + s];
      sdata[threadIdx.x + blockDim.x * 2]
        += sdata[threadIdx.x + blockDim.x * 2 + s];
      sdata[threadIdx.x + blockDim.x * 3]
        += sdata[threadIdx.x + blockDim.x * 3 + s];
      sdata[threadIdx.x + blockDim.x * 4]
        += sdata[threadIdx.x + blockDim.x * 4 + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    ds[blockIdx.x + blockIdx.y * gridDim.x] = sdata[0];
    ds[blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y] =
      sdata[blockDim.x];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 2] =
      sdata[blockDim.x * 2];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 3] =
      sdata[blockDim.x * 3];
    ds[blockIdx.x + blockIdx.y * gridDim.x + (gridDim.x * gridDim.y) * 4] =
      sdata[blockDim.x * 4];
  }
}

// theta31/theta32, which used to finish the block-sum reduction (ds ->
// ds2 -> ds) fully on the device, were removed: AN2::cal_theta() has
// finished that last step with thrust::reduce() for a while and neither
// kernel was called from anywhere, so they were dead code kept around
// from an earlier version of this reduction.
