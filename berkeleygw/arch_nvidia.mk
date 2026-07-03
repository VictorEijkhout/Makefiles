# arch.mk for BerkeleyGW codes
#
# Suitable for Ubuntu 10.10 in parallel.
# Install packages: liblapack-dev, fftw-dev, gfortran, g++, mpi-default-dev
# Build BLACS according to http://www.open-mpi.org/faq/?category=mpi-apps#blacs
# Build ScaLAPACK according to http://www.open-mpi.org/faq/?category=mpi-apps#scalapack,
#
# D. Strubbe
# January 2011, UCB

# VLE just guessing
COMPFLAG  = -DNVIDIA
PARAFLAG  = -DMPI -DOMP

MATHFLAG  =  -DUSESCALAPACK -DUSEFFTW3 -DHDF5

# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

#########################################################################
#   NOTE: This arch.mk is used by a buildslave. The compiler flags are  #
#   optimized to debug the code and not for code performance.           #
#########################################################################

# VLE cpp is found in TACC_GCC_BIN or /usr/bin
# there doesn't seem to be much difference.
# Current options due to Albert
FCPP    = cpp -C -nostdinc

# VLE there are MPI type errors;
# leave out "-Wall -pedantic-errors"
F90free = mpif90 -Mfree

LINK    = mpif90 
FOPTS   = -O3 
FNOOPTS = -O2 
# VLE note the space at the end of this line!
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP  = mpicxx -Wall -std=c++0x
C_COMP   = mpicc  -Wall -std=c99
C_LINK   = mpicxx
C_OPTS   = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
MKLPATH      = ${TACC_NVPL_LIB}

FFTWLIB      =	-Wl,--start-group \
		$(MKLPATH)/libnvpl_intel_lp64.a \
		$(MKLPATH)/libnvpl_intel_thread.a \
		$(MKLPATH)/libnvpl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
FFTWINCLUDE  = ${TACC_NVPL_INC}


LAPACKLIB    = -Wl,--start-group \
		$(MKLPATH)/libnvpl_intel_lp64.a \
		$(MKLPATH)/libnvpl_intel_thread.a \
		$(MKLPATH)/libnvpl_core.a \
		$(MKLPATH)/libnvpl_blacs_intelmpi_lp64.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
SCALAPACKLIB = $(MKLPATH)/libnvpl_scalapack_lp64.a

# HDF5PATH     = /home1/apps/intel25/impi21/phdf5/1.14.6/lib
HDF5PATH = ${TACC_PHDF5_LIB}
HDF5LIB      =	$(HDF5PATH)/libhdf5_hl_fortran.a \
		$(HDF5PATH)/libhdf5_hl.a \
		$(HDF5PATH)/libhdf5_fortran.a \
                $(HDF5PATH)/libhdf5_f90cstub.a \
		$(HDF5PATH)/libhdf5.a \
		/usr/lib64/libsz.so \
		-lz
HDF5INCLUDE  = ${HDF5PATH}/../include
