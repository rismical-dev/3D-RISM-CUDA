#include "rism3d.h"
#include "cuda_check.h"

void RISM3D :: set_cuda () {
  b.x = ce -> grid[0];
  g.x = ce -> grid[1];
  g.y = ce -> grid[2];

  // Every kernel launched with RISM3D's own g/b (cal_rmdft.cu's sum/sum2/
  // cal_dwork/..., and the other cal_*.cu/output_*.cu kernels) as well as
  // AN2's (anderson.cu, which is initialized later from the same ce->grid
  // values in iterate()'s "ma -> initialize(ce, sv)") uses b.x == grid[0]
  // threads per block. grid[0] is chosen for numerical reasons (grid
  // spacing vs. the box size, see cell.cc's MAX_DR check), not to fit
  // inside a CUDA block, so a fine enough grid can quietly ask for more
  // threads per block than the GPU allows. Catching that here, right after
  // read_input() and before any kernel launch in this run, is far easier
  // to debug than the launch failure (or silent no-op, depending on the
  // driver) that would otherwise surface deep inside cal_potential(),
  // iterate(), or output(). This only fires once grid[0] actually exceeds
  // the device's limit, so it changes nothing for any grid size that
  // already works. Because this check runs before AN2::initialize() ever
  // does (set_cuda() is called from RISM3D::initialize(), well before
  // iterate()'s "ma -> initialize(ce, sv)"), and both checks test the same
  // ce->grid[0] against the same device limit, AN2::initialize()'s own
  // copy of this check is now unreachable dead code and was removed.
  int dev;
  AN_CUDA_CHECK(cudaGetDevice(&dev));
  cudaDeviceProp prop;
  AN_CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
  if ((int)b.x > prop.maxThreadsPerBlock) {
    fprintf(stderr,
            "RISM3D::set_cuda: grid[0] = %u exceeds this device's "
            "maxThreadsPerBlock = %d; reduce the grid size or restructure "
            "the launch to split the x-axis across multiple blocks.\n",
            b.x, prop.maxThreadsPerBlock);
    exit(1);
  }
}
