#ifndef AN_CUFFT_CHECK_H
#define AN_CUFFT_CHECK_H

#include <cstdio>
#include <cstdlib>
#include <cufft.h>

// cuFFT has no built-in cufftGetErrorString(), so the mapping is done here.
// (This also matches whatever cuFFT compatibility layer SCALE provides for
// the AMD build, since the cufftResult codes themselves are just an enum.)
inline const char * anCufftErrorString (cufftResult err) {
  switch (err) {
    case CUFFT_SUCCESS:              return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:         return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:         return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:         return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:       return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:          return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:         return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:         return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:       return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_INVALID_DEVICE:       return "CUFFT_INVALID_DEVICE";
    case CUFFT_PARSE_ERROR:          return "CUFFT_PARSE_ERROR";
    case CUFFT_NO_WORKSPACE:         return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:      return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_NOT_SUPPORTED:        return "CUFFT_NOT_SUPPORTED";
    default:                         return "UNKNOWN_CUFFT_ERROR";
  }
}

// Checks the result of a cuFFT call and aborts immediately with the
// file/line and the error name, instead of letting a bad plan or a
// failed exec silently produce garbage (or, on SCALE/AMD, a later,
// unrelated-looking crash).
#define AN_CUFFT_CHECK(call) do {                                           \
    cufftResult an_cufft_err__ = (call);                                    \
    if (an_cufft_err__ != CUFFT_SUCCESS) {                                  \
      fprintf(stderr, "cuFFT error at %s:%d: %s (%d)\n", __FILE__, __LINE__,\
              anCufftErrorString(an_cufft_err__), (int)an_cufft_err__);     \
      exit(1);                                                              \
    }                                                                       \
  } while (0)

#endif // AN_CUFFT_CHECK_H
