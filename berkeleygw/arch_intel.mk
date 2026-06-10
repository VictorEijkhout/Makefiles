# arch.mk for BerkeleyGW codes
#
# Suitable for Ubuntu 10.10 in parallel.
# Install packages: liblapack-dev, fftw-dev, gfortran, g++, mpi-default-dev
# Build BLACS according to http://www.open-mpi.org/faq/?category=mpi-apps#blacs
# Build ScaLAPACK according to http://www.open-mpi.org/faq/?category=mpi-apps#scalapack,
#
# D. Strubbe
# January 2011, UCB

COMPFLAG  = -DINTEL 
PARAFLAG  = -DMPI -DOMP

## VLE we will get to this next: -DUSESCALAPACK
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
F90free = mpif90 -free -qopenmp -no-ipo

LINK    = mpif90 -qopenmp -no-ipo
FOPTS   = -O3 -fp-model source
FNOOPTS = -O2 -fp-model source -no-ipo
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
MKLPATH      = $(MKLROOT)/lib/intel64

FFTWLIB      =	-Wl,--start-group \
		$(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a \
		$(MKLPATH)/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
FFTWINCLUDE  = $(MKLROOT)/include/fftw


LAPACKLIB    = -Wl,--start-group \
		$(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a \
		$(MKLPATH)/libmkl_core.a \
		$(MKLPATH)/libmkl_blacs_intelmpi_lp64.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a

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
