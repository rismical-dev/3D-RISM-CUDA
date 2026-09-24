# Overridable from the command line, e.g.
#   make CUDA_DIR=/opt/cuda-12.4 GPU_ARCH=sm_80
# GPU_ARCH is the -gencode target (sm_90 = H100). Named GPU_ARCH rather
# than ARCH because some clusters' module environments export ARCH
# (e.g. ARCH=x86_64), which ?= would silently pick up.
# For the SCALE/AMD build, point CUDA_DIR at the SCALE target directory
# and set GPU_ARCH to whatever sm_XX that SCALE install maps to gfx942.
CUDA_DIR ?= /usr/local/cuda
GPU_ARCH ?= sm_90

CC 	= g++

CFLAGS 	= -O3 -fomit-frame-pointer -fopenmp -std=c++14

NVCC	= nvcc

NVFLAGS = -O3 -m64 -gencode=arch=$(subst sm_,compute_,$(GPU_ARCH)),code=$(GPU_ARCH) -I$(CUDA_DIR)/include

LDFLAGS = -L$(CUDA_DIR)/lib64
# Runtime API only (cuda* calls) -- no driver API (cu* calls) is used, so
# -lcuda isn't needed; libcudart loads the driver itself at run time.
LIBS = -lcufft -lcudart

OBJS   	= main.o\
	  add_tuv.o\
	  alloc2D.o\
	  alloc3D.o\
	  anderson.o\
	  anderson_cuda.o\
	  cal_ad1.o\
	  cal_ad2.o\
	  cal_Coulomb.o\
	  cal_euv.o\
	  cal_exchem.o\
	  cal_rmdft.o\
	  cal_grad.o\
	  cal_LJ.o\
	  cal_pmv.o\
	  cal_potential.o\
	  cal_pressure.o\
	  cal_qv.o\
	  cal_rms.o\
	  calculate.o\
	  calculate_cuda.o\
	  calloc2D.o\
	  cell.o\
	  dealloc2D.o\
	  dealloc3D.o\
	  excp_fd.o\
	  fft3d.o\
	  initialize.o\
	  initialize_g.o\
	  initialize_tuv.o\
	  iterate.o\
	  output.o\
	  output_ad.o\
	  output_cuv.o\
	  output_euv.o\
	  output_grad.o\
	  output_guv.o\
	  output_huv.o\
	  output_qv.o\
	  output_xmu.o\
	  read_input.o\
	  read_tuv.o\
	  set_ad.o\
	  set_cuda.o\
	  set_fname.o\
	  set_solvent.o\
	  solute.o\
	  solvent.o\
	  solvent_read.o\
	  solvent_spline.o\
	  solvent_spline2.o\
	  spline.o\
	  splint.o\
	  write_tuv.o\

PROGRAM = 3drism-cuda

all:	$(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

# Header dependency tracking: each compile also writes a .d file listing
# every non-system header the source includes (-MMD), plus an empty rule
# per header (-MP) so that deleting/renaming a header doesn't break the
# build. Those .d files are pulled in below, so touching e.g. rism3d.h or
# grid_constants.h now rebuilds every object that includes it. Before
# this, the %.o rules only looked at the .cc/.cu file's own timestamp, so
# a header edit could silently leave stale objects in the link.
DEPFLAGS = -MMD -MP -MF $(@:.o=.d)

%.o: %.cc
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVFLAGS) $(DEPFLAGS) -c $< -o $@

.PHONY: all clean

clean:
	@rm -f $(OBJS) $(OBJS:.o=.d) core

-include $(OBJS:.o=.d)
