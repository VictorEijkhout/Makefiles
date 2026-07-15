# arch.mk for BerkeleyGW codes
#
# Suitable for Ubuntu 10.10 in parallel.
# Install packages: liblapack-dev, fftw-dev, gfortran, g++, mpi-default-dev
# Build BLACS according to http://www.open-mpi.org/faq/?category=mpi-apps#blacs
# Build ScaLAPACK according to http://www.open-mpi.org/faq/?category=mpi-apps#scalapack,
#
# D. Strubbe
# January 2011, UCB

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI
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
# Setting "-x none" leads to a "gfortran -E" error. Hm.
# Leave out the "-C" option fixes it???
FCPP    = cpp 

#F90free = /usr/bin/mpif90 -ffree-form -ffree-line-length-none -Wall
#F90free = mpif90 -ffree-form -ffree-line-length-none -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

# VLE there are MPI type errors;
# leave out "-Wall -pedantic-errors"
F90free = mpif90 -ffree-form -ffree-line-length-none -fbounds-check -std=gnu

# denormal not supported, according to runtime error
# precision is flagged by cputime: http://gcc.gnu.org/onlinedocs/gcc-4.5.3/gfortran/Debugging-Options.html
# underflow is flagged by dlamch from lapack
LINK    = mpif90
FOPTS   = -O3
#FOPTS   = -g -C
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

C_PARAFLAG  = -DPARA
CC_COMP  = mpicxx -Wall -pedantic-errors -std=c++0x
C_COMP  = mpicc -Wall -pedantic-errors -std=c99
C_LINK  = mpicxx
C_OPTS  = -O3
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L${TACC_FFTW3_LIB} -lfftw3
FFTWINCLUDE  = ${TACC_FFTW3_INC}
LAPACKLIB    = -L/usr/lib/ -llapack
BLACSDIR     = /usr/lib
BLACS        = $(BLACSDIR)/blacs_MPI-LINUX-1.a $(BLACSDIR)/blacsF77init_MPI-LINUX-1.a $(BLACSDIR)/blacs_MPI-LINUX-1.a
SCALAPACKLIB = /usr/lib/libscalapack-1.a $(BLACS)

#Ubuntu BLACS/SCALAPACK packages don't work properly sometimes
#BLACS        = -lblacs-openmpi -lblacsF77init-openmpi -lblacs-openmpi
#SCALAPACK    = -lscalapack-openmpi $(BLACS)

#MKL ScaLAPACK fails completely in complex version
#MKLPATH      = /auto/opt/intel/mkl/lib/intel64
#LAPACKLIB    = -Wl,--start-group $(MKLPATH)/libmkl_gf_lp64.a $(MKLPATH)/libmkl_sequential.a \
               $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread
#SCALAPACKLIB = $(MKLPATH)/libmkl_scalapack_lp64.a 

#need to export MPIEXEC=/usr/bin/mpirun if this is not default in `which mpiexec`
TESTSCRIPT = make check-parallel

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
