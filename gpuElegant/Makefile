TOP=../../..
include $(TOP)/configure/CONFIG
include $(TOP)/configure/RULES

NVCC = $(wildcard /usr/local/cuda/bin/nvcc /usr/local/cuda-8.0/bin/nvcc /usr/local/cuda-7.5/bin/nvcc)

MPI_INCLUDE = $(wildcard /usr/include/mpich-x86_64 )

ifeq ($(MPI),1)
all: elegantto_gpu.o elegant_gpu.o Pelegant_gpu.o
else
all: elegantto_gpu.o elegant_gpu.o
endif

clean::
	rm -f *.o

OBJNAMES = gpu_base.o \
	gpu_bin_time_distribution.o \
	gpu_bin_transverse_distribution.o \
	gpu_compute_centroids.o \
	gpu_convolve_arrays.o \
	gpu_csbend.o \
	gpu_final_props.o \
	gpu_funcs.o \
	gpu_kahan.o \
	gpu_killParticles.o \
	gpu_limit_amplitudes.o \
	gpu_lrwake.o \
	gpu_lsc.o \
	gpu_matrix.o \
	gpu_matter.o \
	gpu_memory_management.o \
	gpu_multipole.o \
	gpu_reductions.o \
	gpu_simple_rfca.o \
	gpu_trwake.o \
	gpu_wake.o \
	printCudaDeviceProp.o

MPIOBJNAMES = mpi_gpu_base.o \
	mpi_gpu_bin_time_distribution.o \
	mpi_gpu_bin_transverse_distribution.o \
	mpi_gpu_compute_centroids.o \
	mpi_gpu_convolve_arrays.o \
	mpi_gpu_csbend.o \
	mpi_gpu_final_props.o \
	mpi_gpu_funcs.o \
	mpi_gpu_kahan.o \
	mpi_gpu_killParticles.o \
	mpi_gpu_limit_amplitudes.o \
	mpi_gpu_lrwake.o \
	mpi_gpu_lsc.o \
	mpi_gpu_matrix.o \
	mpi_gpu_matter.o \
	mpi_gpu_memory_management.o \
	mpi_gpu_multipole.o \
	mpi_gpu_reductions.o \
	mpi_gpu_simple_rfca.o \
	mpi_gpu_trwake.o \
	mpi_gpu_wake.o \
	mpi_printCudaDeviceProp.o

ELEGANTTOOBJNAMES = gpu_base.o gpu_matrix.o gpu_memory_management.o

HCUNAMES = gpu_fringe.hcu \
	   gpu_interp.hcu \
	   gpu_kahan.hcu \
	   gpu_killParticles.hcu \
	   gpu_matter.hcu \
	   gpu_particle_accessor.hcu \
	   gpu_particle_reduction.hcu \
	   gpu_particle_template_function.hcu \
	   gpu_reductions.hcu

HNAMES = gpu_base.h \
	 gpu_bin_time_distribution.h \
	 gpu_bin_transverse_distribution.h \
	 gpu_compute_centroids.h \
	 gpu_convolve_arrays.h \
	 gpu_csbend.h \
	 gpu_final_props.h \
	 gpu_funcs.h \
	 gpu_kahan.h \
	 gpu_limit_amplitudes.h \
	 gpu_lrwake.h \
	 gpu_lsc.h \
	 gpu_matrix.h \
	 gpu_matter.h \
	 gpu_memory_management.h \
	 gpu_multipole.h \
	 gpu_reductions.h \
	 gpu_simple_rfca.h \
	 gpu_track.h \
	 gpu_trwake.h \
	 gpu_verify.h \
	 gpu_wake.h \
	 printCudaDeviceProp.h


CFLAGS = -c -ccbin /usr/bin/cc -m64 -DHAVE_CONFIG_H -DzLib -DHAVE_LAPACK -DLAPACK -DUSE_KAHAN -DIEEE_MATH -DMINIMIZE_MEMORY -DWL=32 -DRPN_SUPPORT -DUNIX -Dlinux -D_X86_64_ -DHAVE_GPU=1 -DUSE_GSL -Xcompiler ,\"-ansi\",\"-pipe\",\"-O3\",\"-DNDEBUG\",\"-msse2\" -Xcompiler ,\"-Wno-write-strings\" --fmad=true --Werror cross-execution-space-call -O3 --use_fast_math --generate-code arch=compute_52,code=sm_52 -DNVCC -I. -I.. -I$(EPICS_EXTENSIONS_INCLUDE) -I$(EPICS_EXTENSIONS)/src/SDDS/lzma

#CFLAGS += -DCUDA_DEVICE_COUNT=<count>
#CFLAGS += -DDEBUG -DCUDA_DEBUG -DMPI_DEBUG -DDEBUG_APERTURE -DDEBUG_IGF
#CFLAGS += -DDEBUG_FIDUCIALIZATION


ifeq ($(VERIFY),1)
CFLAGS += -DGPU_VERIFY
OBJNAMES += gpu_verify.o
MPIOBJNAMES += mpi_gpu_verify.o
ELEGANTTOOBJNAMES += gpu_verify.o
endif

MPI_CFLAGS=-DUSE_MPI -DSDDS_MPI_IO -I$(MPI_INCLUDE)

%.o: %.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $<

mpi_gpu_base.o: gpu_base.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_base.o gpu_base.cu

mpi_gpu_bin_time_distribution.o: gpu_bin_time_distribution.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_bin_time_distribution.o gpu_bin_time_distribution.cu

mpi_gpu_bin_transverse_distribution.o: gpu_bin_transverse_distribution.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_bin_transverse_distribution.o gpu_bin_transverse_distribution.cu

mpi_gpu_compute_centroids.o: gpu_compute_centroids.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_compute_centroids.o gpu_compute_centroids.cu

mpi_gpu_convolve_arrays.o: gpu_convolve_arrays.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_convolve_arrays.o gpu_convolve_arrays.cu

mpi_gpu_csbend.o: gpu_csbend.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_csbend.o gpu_csbend.cu

mpi_gpu_final_props.o: gpu_final_props.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_final_props.o gpu_final_props.cu

mpi_gpu_funcs.o: gpu_funcs.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_funcs.o gpu_funcs.cu

mpi_gpu_kahan.o: gpu_kahan.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_kahan.o gpu_kahan.cu

mpi_gpu_killParticles.o: gpu_killParticles.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_killParticles.o gpu_killParticles.cu

mpi_gpu_limit_amplitudes.o: gpu_limit_amplitudes.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_limit_amplitudes.o gpu_limit_amplitudes.cu

mpi_gpu_lrwake.o: gpu_lrwake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_lrwake.o gpu_lrwake.cu

mpi_gpu_lsc.o: gpu_lsc.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_lsc.o gpu_lsc.cu

mpi_gpu_matrix.o: gpu_matrix.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_matrix.o gpu_matrix.cu

mpi_gpu_matter.o: gpu_matter.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_matter.o gpu_matter.cu

mpi_gpu_memory_management.o: gpu_memory_management.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_memory_management.o gpu_memory_management.cu

mpi_gpu_multipole.o: gpu_multipole.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_multipole.o gpu_multipole.cu

mpi_gpu_reductions.o: gpu_reductions.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_reductions.o gpu_reductions.cu

mpi_gpu_simple_rfca.o: gpu_simple_rfca.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_simple_rfca.o gpu_simple_rfca.cu

mpi_gpu_trwake.o: gpu_trwake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_trwake.o gpu_trwake.cu

mpi_gpu_wake.o: gpu_wake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_wake.o gpu_wake.cu

mpi_printCudaDeviceProp.o: printCudaDeviceProp.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_printCudaDeviceProp.o printCudaDeviceProp.cu

mpi_gpu_verify.o: gpu_verify.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_verify.o gpu_verify.cu


elegantto_gpu.o: $(ELEGANTTOOBJNAMES)
	$(NVCC) --gpu-architecture=sm_52 --device-link \
	        $(ELEGANTTOOBJNAMES) --output-file elegantto_gpu.o -lcudadevrt -lcudart

elegant_gpu.o: $(OBJNAMES)
	$(NVCC) --generate-code arch=compute_52,code=sm_52 --device-link \
	        $(OBJNAMES) --output-file elegant_gpu.o -lcudadevrt -lcudart

Pelegant_gpu.o: $(MPIOBJNAMES)
	$(NVCC) --gpu-architecture=sm_52 --device-link \
	        $(MPIOBJNAMES) --output-file Pelegant_gpu.o -lcudadevrt -lcudart


