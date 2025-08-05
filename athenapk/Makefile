################################################################
####
#### Makefile for AthenaK (not Athenapk) installation
####
################################################################

PACKAGE = athenak
URL = https://github.com/IAS-Astrophysics/athenak
# parthenon-hpc-lab/athenapk
ABOUT = Performance-portable version of the Athena++ astrophysical AMR-MHD code using Kokkos
PACKAGEVERSION = git

MODE = mpi
MODULES = phdf5 cuda

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## cmake version
##
CPPSTANDARD = 20
CUDAFLAGS = -allow-unsupported-compiler -Xcompiler -std=gnu++${CPPSTANDARD}
.PHONY: outer_config
outer_config :
	export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure 

CMAKEFLAGS = \
    -D Athena_ENABLE_MPI=ON \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_CUDA=On \
    -D CMAKE_CUDA_FLAGS='${CUDAFLAGS}' \
    -D Kokkos_ARCH_HOPPER90=ON \
    -D CMAKE_CXX_STANDARD=${CPPSTANDARD}
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : outer_config build"
	@echo "Do Not Use 'config' target"
default_install : outer_config build

GITREPO = https://github.com/IAS-Astrophysics/athenak.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

