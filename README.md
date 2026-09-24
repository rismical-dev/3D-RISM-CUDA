### INSTALL

Requirements: CUDA Toolkit (nvcc, cuFFT) and g++ with OpenMP support.
The AMD build uses the SCALE compiler in place of the CUDA Toolkit.

Build with default settings (CUDA in /usr/local/cuda, target sm_90 = H100):

<pre>
make
</pre>

The CUDA location and GPU target can be set on the command line
instead of editing the Makefile:

<pre>
make CUDA_DIR=/opt/cuda-12.4 GPU_ARCH=sm_80
</pre>

<pre>
Variable   Default           Meaning
--------------------------------------------------------------------
CUDA_DIR   /usr/local/cuda   CUDA (or SCALE) installation directory
GPU_ARCH   sm_90             GPU target passed to nvcc -gencode
</pre>

Both use `?=`, so an environment variable of the same name also takes
effect. If the build picks up an unexpected path, check `echo $CUDA_DIR`.

For AMD GPUs with SCALE, set `CUDA_DIR` to the SCALE target directory and
`GPU_ARCH` to the `sm_XX` that your SCALE installation maps to the AMD GPU
(e.g. gfx942 for MI300).

Header dependencies are tracked automatically (`.d` files), so after
editing a header only the files that include it are recompiled.
Run `make clean` to remove objects and dependency files; do this once
after updating from an older Makefile.

The program uses only the CUDA runtime API, so it links against
`libcudart` and `libcufft`; `libcuda` (the driver library) is not needed,
and the program can be built on nodes without a GPU driver.

### PERFORMANCE

Input: test4.inp

Solute: ProteinG (846 atoms)

Solvent: Water (Tip3p model)

Cell: 128^3 Angstrom^3 (256^3 grids)

<pre>
GPU      Time [s]
-----------------
K20C     43.8
P100     14.5
Titan V  12.5
V100      9.9
</pre>

### REFERENCE
1. Y. Maruyama, and N. Yoshida, "RISMiCal: A software package to perform fast RISM/3D-RISM calculations," J. Comput. Chem., (2024) 45, 1470-1482 (DOI: 10.1002/jcc.27340)

2. Y. Maruyama and F. Hirata, "Modified Anderson method for accelerating 3D-RISM calculations using graphics processing unit," J. Chem. Theory Comput., (2012) 8, 3015-3021 (DOI: 10.1021/ct300355r)

3. Y. Maruyama, and N. Yoshida, "High-Precision Solvation Free Energy Calculation via Multi-Input Linear Correction in 3D-RISM Theory",  J. Chem. Theory Comput., (2026) 22, 3596–3612(DOI: 10.1021/acs.jctc.5c02013)