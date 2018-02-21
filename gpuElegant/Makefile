TOP=../../..
include $(TOP)/configure/CONFIG

ifeq ($(EPICS_HOST_ARCH),windows-x64)
WIN32 = 1
endif

ifdef WIN32
NVCC = nvcc
OBJ = obj
else
NVCC = $(firstword $(wildcard /usr/bin/nvcc /usr/local/cuda-9.1/bin/nvcc /usr/local/cuda-9.0/bin/nvcc /usr/local/cuda-8.0/bin/nvcc /usr/local/cuda-7.5/bin/nvcc /usr/local/cuda/bin/nvcc /Developer/NVIDIA/CUDA-9.1/bin/nvcc /Developer/NVIDIA/CUDA-9.0/bin/nvcc /Developer/NVIDIA/CUDA-8.0/bin/nvcc /soft/cudatoolkit/6.5.14/bin/nvcc))
OBJ = o
endif

MPI_INCLUDE = $(wildcard /usr/include/mpich-3.2-x86_64 /software/mvapich2-intel-psm-1.9.5/include )

CUDA_HARDWARE_52=52
CUDA_HARDWARE_60=60
#Check if we are on Blues
ifeq (blogin, $(findstring blogin, $(HOSTNAME)))
CUDA_HARDWARE=35
endif

OBJNAMES = gpu_base.$(OBJ) \
	gpu_bin_time_distribution.$(OBJ) \
	gpu_bin_transverse_distribution.$(OBJ) \
	gpu_compute_centroids.$(OBJ) \
	gpu_convolve_arrays.$(OBJ) \
	gpu_csbend.$(OBJ) \
	gpu_final_props.$(OBJ) \
	gpu_funcs.$(OBJ) \
	gpu_kahan.$(OBJ) \
	gpu_killParticles.$(OBJ) \
	gpu_limit_amplitudes.$(OBJ) \
	gpu_lrwake.$(OBJ) \
	gpu_lsc.$(OBJ) \
	gpu_matrix.$(OBJ) \
	gpu_matter.$(OBJ) \
	gpu_memory_management.$(OBJ) \
	gpu_multipole.$(OBJ) \
	gpu_reductions.$(OBJ) \
	gpu_simple_rfca.$(OBJ) \
	gpu_trwake.$(OBJ) \
	gpu_wake.$(OBJ) \
	printCudaDeviceProp.$(OBJ)

MPIOBJNAMES = mpi_gpu_base.$(OBJ) \
	mpi_gpu_bin_time_distribution.$(OBJ) \
	mpi_gpu_bin_transverse_distribution.$(OBJ) \
	mpi_gpu_compute_centroids.$(OBJ) \
	mpi_gpu_convolve_arrays.$(OBJ) \
	mpi_gpu_csbend.$(OBJ) \
	mpi_gpu_final_props.$(OBJ) \
	mpi_gpu_funcs.$(OBJ) \
	mpi_gpu_kahan.$(OBJ) \
	mpi_gpu_killParticles.$(OBJ) \
	mpi_gpu_limit_amplitudes.$(OBJ) \
	mpi_gpu_lrwake.$(OBJ) \
	mpi_gpu_lsc.$(OBJ) \
	mpi_gpu_matrix.$(OBJ) \
	mpi_gpu_matter.$(OBJ) \
	mpi_gpu_memory_management.$(OBJ) \
	mpi_gpu_multipole.$(OBJ) \
	mpi_gpu_reductions.$(OBJ) \
	mpi_gpu_simple_rfca.$(OBJ) \
	mpi_gpu_trwake.$(OBJ) \
	mpi_gpu_wake.$(OBJ) \
	mpi_printCudaDeviceProp.$(OBJ)

ELEGANTTOOBJNAMES = gpu_base.$(OBJ) gpu_matrix.$(OBJ) gpu_memory_management.$(OBJ)

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


ifdef WIN32
CFLAGS = -c -ccbin cl -m64 -DHAVE_CONFIG_H -DzLib -DHAVE_LAPACK -DLAPACK -DUSE_KAHAN -DIEEE_MATH -DMINIMIZE_MEMORY -DWL=32 -DRPN_SUPPORT -DUNIX -Dlinux -D_X86_64_ -DHAVE_GPU=1 -DUSE_GSL -D_FORCE_INLINES -Xcompiler ,\"-Ox\",\"-DNDEBUG\",\"-MD\" --fmad=true --Werror cross-execution-space-call -O3 --use_fast_math --generate-code arch=compute_52,code=sm_52 --generate-code arch=compute_60,code=sm_60 -DNVCC -I. -I.. -I$(EPICS_EXTENSIONS_INCLUDE) -I$(EPICS_EXTENSIONS)/src/SDDS/lzma
else
CFLAGS = -c -ccbin /usr/bin/cc -m64 -DHAVE_CONFIG_H -DzLib -DHAVE_LAPACK -DLAPACK -DUSE_KAHAN -DIEEE_MATH -DMINIMIZE_MEMORY -DWL=32 -DRPN_SUPPORT -DUNIX -Dlinux -D_X86_64_ -DHAVE_GPU=1 -DUSE_GSL -D_FORCE_INLINES -Xcompiler ,\"-ansi\",\"-pipe\",\"-O3\",\"-DNDEBUG\",\"-msse2\" -Xcompiler ,\"-Wno-write-strings\" --fmad=true --Werror cross-execution-space-call -O3 --use_fast_math --generate-code arch=compute_$(CUDA_HARDWARE),code=sm_$(CUDA_HARDWARE) -DNVCC -I. -I.. -I$(EPICS_EXTENSIONS_INCLUDE) -I$(EPICS_EXTENSIONS)/src/SDDS/lzma
CFLAGS = -c -ccbin /usr/bin/cc -m64 -DHAVE_CONFIG_H -DzLib -DHAVE_LAPACK -DLAPACK -DUSE_KAHAN -DIEEE_MATH -DMINIMIZE_MEMORY -DWL=32 -DRPN_SUPPORT -DUNIX -Dlinux -D_X86_64_ -DHAVE_GPU=1 -DUSE_GSL -D_FORCE_INLINES -Xcompiler ,\"-ansi\",\"-pipe\",\"-O3\",\"-DNDEBUG\",\"-msse2\" -Xcompiler ,\"-Wno-write-strings\" --fmad=true --Werror cross-execution-space-call -O3 --use_fast_math --generate-code arch=compute_$(CUDA_HARDWARE_52),code=sm_$(CUDA_HARDWARE_52) --generate-code arch=compute_$(CUDA_HARDWARE_60),code=sm_$(CUDA_HARDWARE_60) -DNVCC -I. -I.. -I$(EPICS_EXTENSIONS_INCLUDE) -I$(EPICS_EXTENSIONS)/src/SDDS/lzma
endif

#CFLAGS += -DCUDA_DEVICE_COUNT=<count>
#CFLAGS += -DDEBUG -DCUDA_DEBUG -DMPI_DEBUG -DDEBUG_APERTURE -DDEBUG_IGF
#CFLAGS += -DDEBUG_FIDUCIALIZATION


ifeq ($(VERIFY),1)
CFLAGS += -DGPU_VERIFY
OBJNAMES += gpu_verify.$(OBJ)
MPIOBJNAMES += mpi_gpu_verify.$(OBJ)
ELEGANTTOOBJNAMES += gpu_verify.$(OBJ)
endif

MPI_CFLAGS=-DUSE_MPI -DSDDS_MPI_IO -I$(MPI_INCLUDE)


include $(TOP)/configure/RULES

ifeq ($(MPI),1)
all: elegantto_gpu.$(OBJ) elegant_gpu.$(OBJ) Pelegant_gpu.$(OBJ)
else
all: elegantto_gpu.$(OBJ) elegant_gpu.$(OBJ)
endif

ifdef BASE_3_15
clean:
else
clean::
endif
	rm -f *.$(OBJ)


%.$(OBJ): %.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $<

mpi_gpu_base.$(OBJ): gpu_base.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_base.$(OBJ) gpu_base.cu

mpi_gpu_bin_time_distribution.$(OBJ): gpu_bin_time_distribution.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_bin_time_distribution.$(OBJ) gpu_bin_time_distribution.cu

mpi_gpu_bin_transverse_distribution.$(OBJ): gpu_bin_transverse_distribution.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_bin_transverse_distribution.$(OBJ) gpu_bin_transverse_distribution.cu

mpi_gpu_compute_centroids.$(OBJ): gpu_compute_centroids.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_compute_centroids.$(OBJ) gpu_compute_centroids.cu

mpi_gpu_convolve_arrays.$(OBJ): gpu_convolve_arrays.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_convolve_arrays.$(OBJ) gpu_convolve_arrays.cu

mpi_gpu_csbend.$(OBJ): gpu_csbend.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_csbend.$(OBJ) gpu_csbend.cu

mpi_gpu_final_props.$(OBJ): gpu_final_props.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_final_props.$(OBJ) gpu_final_props.cu

mpi_gpu_funcs.$(OBJ): gpu_funcs.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_funcs.$(OBJ) gpu_funcs.cu

mpi_gpu_kahan.$(OBJ): gpu_kahan.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_kahan.$(OBJ) gpu_kahan.cu

mpi_gpu_killParticles.$(OBJ): gpu_killParticles.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_killParticles.$(OBJ) gpu_killParticles.cu

mpi_gpu_limit_amplitudes.$(OBJ): gpu_limit_amplitudes.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_limit_amplitudes.$(OBJ) gpu_limit_amplitudes.cu

mpi_gpu_lrwake.$(OBJ): gpu_lrwake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_lrwake.$(OBJ) gpu_lrwake.cu

mpi_gpu_lsc.$(OBJ): gpu_lsc.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_lsc.$(OBJ) gpu_lsc.cu

mpi_gpu_matrix.$(OBJ): gpu_matrix.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_matrix.$(OBJ) gpu_matrix.cu

mpi_gpu_matter.$(OBJ): gpu_matter.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_matter.$(OBJ) gpu_matter.cu

mpi_gpu_memory_management.$(OBJ): gpu_memory_management.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_memory_management.$(OBJ) gpu_memory_management.cu

mpi_gpu_multipole.$(OBJ): gpu_multipole.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_multipole.$(OBJ) gpu_multipole.cu

mpi_gpu_reductions.$(OBJ): gpu_reductions.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_reductions.$(OBJ) gpu_reductions.cu

mpi_gpu_simple_rfca.$(OBJ): gpu_simple_rfca.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_simple_rfca.$(OBJ) gpu_simple_rfca.cu

mpi_gpu_trwake.$(OBJ): gpu_trwake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_trwake.$(OBJ) gpu_trwake.cu

mpi_gpu_wake.$(OBJ): gpu_wake.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_wake.$(OBJ) gpu_wake.cu

mpi_printCudaDeviceProp.$(OBJ): printCudaDeviceProp.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_printCudaDeviceProp.$(OBJ) printCudaDeviceProp.cu

mpi_gpu_verify.$(OBJ): gpu_verify.cu $(HCUNAMES) $(HNAMES)
	$(NVCC) $(CFLAGS) $(MPI_CFLAGS) --output-file mpi_gpu_verify.$(OBJ) gpu_verify.cu


elegantto_gpu.$(OBJ): $(ELEGANTTOOBJNAMES)
	$(NVCC) --generate-code arch=compute_$(CUDA_HARDWARE_52),code=sm_$(CUDA_HARDWARE_52) --generate-code arch=compute_$(CUDA_HARDWARE_60),code=sm_$(CUDA_HARDWARE_60) --device-link \
	        $(ELEGANTTOOBJNAMES) --output-file elegantto_gpu.$(OBJ) -lcudadevrt -lcudart

elegant_gpu.$(OBJ): $(OBJNAMES)
	$(NVCC) --generate-code arch=compute_$(CUDA_HARDWARE_52),code=sm_$(CUDA_HARDWARE_52) --generate-code arch=compute_$(CUDA_HARDWARE_60),code=sm_$(CUDA_HARDWARE_60) --device-link \
	        $(OBJNAMES) --output-file elegant_gpu.$(OBJ) -lcudadevrt -lcudart

Pelegant_gpu.$(OBJ): $(MPIOBJNAMES)
	$(NVCC) --generate-code arch=compute_$(CUDA_HARDWARE_52),code=sm_$(CUDA_HARDWARE_52) --generate-code arch=compute_$(CUDA_HARDWARE_60),code=sm_$(CUDA_HARDWARE_60) --device-link \
	        $(MPIOBJNAMES) --output-file Pelegant_gpu.$(OBJ) -lcudadevrt -lcudart


