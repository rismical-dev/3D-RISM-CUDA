#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"
#include "cuda_check.h"

void RISM3D :: iterate(int cu) {
  void alloc2D (std::vector <double *> &, int, int);
  void calloc2D (std::vector <std::complex <double> *> &, int, int);

  double cf, cuf;

  calloc2D (guv, sv -> natv, ce -> ngrid);
  calloc2D (huv, sv -> natv, ce -> ngrid);
  alloc2D (tuv, sv -> natv, ce -> ngrid);

  AN_CUDA_CHECK(cudaMalloc(&dguv, ce -> ngrid * sv -> natv * sizeof(double2)));
  AN_CUDA_CHECK(cudaMalloc(&dhuv, ce -> ngrid * sv -> natv * sizeof(double2)));
  AN_CUDA_CHECK(cudaMalloc(&dt, ce -> ngrid * sv -> natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&dtr, ce -> ngrid * sv -> natv * sizeof(double)));
  AN_CUDA_CHECK(cudaMalloc(&ds, ce -> grid[1] * ce -> grid[2] * sizeof(double)));

  ma -> initialize (ce, sv);
  fft -> initialize (ce);

  std::ifstream in_file ;
  in_file.open((fname + exttuv).c_str());
  bool saved = in_file.is_open();
  in_file.close();

  if (saved) {
    read_tuv();
    cu = 0;
    cf = 1.0;
  } else if (cu == 0) {
    initialize_tuv();
    cf = 1.0;
  } else {
    AN_CUDA_CHECK(cudaMemset(dt, 0, ce -> ngrid * sv -> natv * sizeof(double)));
    cuf = 1.0 /cu;
    cf = 0.0;
  }
  for (int c = 0; c <= cu; ++c) {
    if (c > 0) {
      cf += cuf;
      if (cf > 1.0) cf = 1.0;
      add_tuv(cuf);
    }
    std::cout << "relaxing 3D UV RISM: Charge Up Factor = " << cf << std::endl;
    bool conver = false;
    bool diverge = false;
    for (int istep = 1; istep <= co -> maxstep; ++istep) {
      calculate(cf);
      double rms = cal_rms ();
      diverge = !isfinite(rms);
      if (diverge) {
	break;
      }
      if (rms <= co -> convergence) {
	conver = true;
      } else {
	ma -> calculate (dt, dtr);
      }
      std::cout << " Step = " << istep << " Reside = " << rms << std::endl;
      if (co -> ksave > 0 && istep % co -> ksave == 0) {
	write_tuv();
      }
      if (conver) {
	if (co -> ksave != 0 && c == cu) {
	  write_tuv();
	}
	break;
      }
    }
    if (diverge) {
      // Unlike "reached the step limit without converging" (which still
      // has finite, merely-incomplete dt/dtr, so it's treated as
      // best-effort and falls through to output()), a diverged rms means
      // dt/dtr are already NaN/Inf. Continuing to the next charge-up step
      // (a higher cf, using this corrupted state as its starting point)
      // cannot recover, and letting the caller reach output() afterward
      // would silently write NaN-filled huv/guv/tuv files. Exiting here,
      // like every other unrecoverable condition in this codebase, makes
      // the failure unambiguous instead.
      std::cout << "Calculation diverged." << std::endl;
      exit(1);
    } else if (!conver) {
      std::cout << "3D UV RISM: reached limit # of relaxation steps: "
	   << co -> maxstep << std::endl;
      break;
    }
  }
  // Synchronous, since huv/guv (host buffers) are read back by output()
  // right after iterate() returns: an async copy here would let that read
  // race the D2H transfer (same class of bug fixed in initialize_g.cu).
  for (int iv = 0; iv < sv -> natv; ++iv) {
    AN_CUDA_CHECK(cudaMemcpy(huv[iv], dhuv + (iv * ce -> ngrid),
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault));
    AN_CUDA_CHECK(cudaMemcpy(guv[iv], dguv + (iv * ce -> ngrid),
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault));
  }
  delete ma;
//  delete fft;
}
