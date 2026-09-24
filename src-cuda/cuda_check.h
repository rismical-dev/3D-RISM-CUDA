#ifndef AN_CUDA_CHECK_H
#define AN_CUDA_CHECK_H

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

// Checks the result of a CUDA/HIP runtime call (or cudaGetLastError() /
// cudaPeekAtLastError() called right after a kernel launch) and aborts
// immediately with the file/line and the real error string, instead of
// letting a silent asynchronous failure surface later as an unrelated-
// looking error inside some other CUDA/HIP/thrust call (which is what
// makes these bugs look like they "happen at output" when the real fault
// is earlier).
#define AN_CUDA_CHECK(call) do {                                            \
    cudaError_t an_cuda_err__ = (call);                                     \
    if (an_cuda_err__ != cudaSuccess) {                                     \
      fprintf(stderr, "CUDA/HIP error at %s:%d: %s\n", __FILE__, __LINE__,  \
              cudaGetErrorString(an_cuda_err__));                           \
      exit(1);                                                              \
    }                                                                       \
  } while (0)

#endif // AN_CUDA_CHECK_H
