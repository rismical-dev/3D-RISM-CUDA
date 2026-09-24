#include "rism3d.h"
#include "cuda_check.h"

void RISM3D :: cal_potential(std::string esp) {
  cal_LJ();
  cal_Coulomb(esp);

  // dgv (built in initialize_g()) is only used by cal_Coulomb.cu's fk
  // kernel, so it's safe to free once cal_LJ()/cal_Coulomb() have run.
  AN_CUDA_CHECK(cudaFree(dgv));
}
