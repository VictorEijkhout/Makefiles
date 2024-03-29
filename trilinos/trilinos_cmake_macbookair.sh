cmake -Wno-dev \
  -D CMAKE_INSTALL_PREFIX:PATH=${installdir} \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D Trilinos_VERBOSE_CONFIGURE=ON \
  -D CMAKE_VERBOSE_MAKEFILE=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  -D Trilinos_ENABLE_Export_Makefiles:BOOL=ON \
  -D Trilinos_ENABLE_Fortran:BOOL=ON \
  \
  -D CMAKE_C_COMPILER=${CC} \
  -D CMAKE_CXX_COMPILER=${CXX} \
  -D CMAKE_Fortran_COMPILER=${FC} \
  -D CMAKE_C_FLAGS:STRING="${COPTFLAGS}" \
  -D CMAKE_CXX_FLAGS:STRING="${COPTFLAGS} -DMPICH_SKIP_MPICXX" \
  -D Trilinos_EXTRA_LINK_FLAGS=${PYTHON_LOAD_FLAG} \
  \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_BASE_DIR=${LMOD_IMPI_DIR}/intel64 \
  -D MPI_INCLUDE_DIRS=${LMOD_IMPI_INC} \
  -D MPI_EXEC:FILEPATH="/opt/apps/xalt/0.6/bin/ibrun" \
  -D TPL_MPI_BASE_DIR=${LMOD_IMPI_DIR}intel64 \
  -D TPL_MPI_INCLUDE_DIRS=${LMOD_IMPI_INC} \
  -D TPL_ENABLE_GLM=OFF \
  -D TPL_ENABLE_Matio=OFF \
  \
  -D TPL_ENABLE_BLAS=ON \
  -D BLAS_INCLUDE_DIRS:PATH="${LMOD_BLIS_INC}" \
  -D BLAS_LIBRARY_DIRS:PATH="${LMOD_BLIS_LIB}" \
  -D BLAS_LIBRARY_NAMES:STRING="blis" \
  -D LAPACK_INCLUDE_DIRS:PATH="${LMOD_BLIS_INC}" \
  -D LAPACK_LIBRARY_DIRS:PATH="${LMOD_BLIS_LIB}" \
  -D LAPACK_LIBRARY_NAMES:STRING="blis" \
  \
  -D TPL_ENABLE_HDF5:BOOL=${HAS_HDF5} \
  -D HDF5_INCLUDE_DIRS:PATH=$LMOD_HDF5_INC    \
  -D HDF5_LIBRARY_DIRS:PATH=$LMOD_HDF5_LIB    \
  -D TPL_ENABLE_Netcdf:BOOL=${HAS_NETCDF} \
  -D Netcdf_INCLUDE_DIRS:PATH=$LMOD_NETCDF_INC    \
  -D Netcdf_LIBRARY_DIRS:PATH=$LMOD_NETCDF_LIB    \
  \
  -D Tpetra_INST_DOUBLE:BOOL=ON \
  -D Tpetra_INST_FLOAT:BOOL=OFF \
  -D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF \
  -D Tpetra_INST_COMPLEX_DOUBLE:BOOL=OFF \
  -D Tpetra_INST_INT_LONG:BOOL=OFF \
  -D Tpetra_INST_INT_UNSIGNED:BOOL=OFF \
  \
  -D TPL_ENABLE_Boost:BOOL=OFF \
  -D Boost_INCLUDE_DIRS:PATH=$LMOD_BOOST_INC      \
  -D Boost_LIBRARY_DIRS:PATH=$LMOD_BOOST_LIB      \
  -D TPL_ENABLE_BoostLib:BOOL=OFF \
  -D BoostLib_INCLUDE_DIRS:PATH=$LMOD_BOOST_INC      \
  -D BoostLib_LIBRARY_DIRS:PATH=$LMOD_BOOST_LIB      \
  \
  -D TPL_ENABLE_MUMPS:BOOL=${HAS_MUMPS} \
  -D MUMPS_INCLUDE_DIRS=${LMOD_MUMPS_INC} \
  -D MUMPS_LIBRARY_DIRS=${LMOD_MUMPS_LIB} \
  -D MUMPS_LIBRARY_NAMES:STRING=${MUMPSLIBNAMES} \
  -D VLE_TPL_MUMPS_LIBRARIES:STRING=${MUMPSLIBS} \
  -D TPL_ParMETIS_LIBRARIES="${PARMETISLIBS}" \
  -D TPL_ParMETIS_INCLUDE_DIRS=${LMOD_PARMETIS_INC} \
  \
  -D TPL_ENABLE_yaml-cpp:BOOL=${HAS_YAMLCPP} \
  -D yaml-cpp_INCLUDE_DIRS:PATH=${LMOD_YAMLCPP_INC} \
  -D yaml-cpp_LIBRARY_DIRS:PATH=${LMOD_YAMLCPP_LIB} \
  \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
      -D Trilinos_ENABLE_Amesos2:BOOL=ON \
      -D Amesos2_ENABLE_KLU2:BOOL=ON \
      -D Amesos2_ENABLE_Basker:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:Bool=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:Bool=ON \
  -D Trilinos_ENABLE_EpetraExt:Bool=ON \
  -D                 Epetra_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_FEI:Bool=ON \
  -D Trilinos_ENABLE_Ifpack:Bool=ON \
      -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
  -D Trilinos_ENABLE_Intrepid:BOOL=ON \
      -D Trilinos_ENABLE_Intrepid2:BOOL=ON \
      -D Intrepid_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_Isorropia:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
      -D ML_TAKES_SUPERLU_LESS_THAN_5=TRUE \
      -D ML_ENABLE_SuperLU:BOOL=OFF \
  -D Trilinos_ENABLE_MOOCHO:BOOL=ON \
  -D Trilinos_ENABLE_MueLu:BOOL=${HAS_MUELU} \
  -D                 MueLu_ENABLE_Tutorial:BOOL=OFF \
  -D                 MueLu_ENABLE_EXAMPLES:BOOL=OFF \
  -D Trilinos_ENABLE_NOX=ON \
  -D                 NOX_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Pamgen:Bool=ON \
  -D Trilinos_ENABLE_Panzer:Bool=ON \
  -D Trilinos_ENABLE_Phalanx:BOOL=ON \
      -D Phalanx_EXPLICIT_TEMPLATE_INSTANTIATION=ON \
      -D Phalanx_ENABLE_EXAMPLES=OFF \
  -D Trilinos_ENABLE_Piro:BOOL=ON \
  -D Trilinos_ENABLE_Rythmos:BOOL=ON \
  -D Trilinos_ENABLE_Sacado:Bool=ON \
  -D Trilinos_ENABLE_SEACAS:BOOL=${HAS_SEACAS} \
      -D Trilinos_ENABLE_SEACASIoss:BOOL=${HAS_SEACAS} \
      -D Trilinos_ENABLE_SEACASBlot:BOOL=${HAS_SEACAS} \
      -D Trilinos_ENABLE_SEACASExodus:BOOL=${HAS_SEACAS} \
  -D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON \
  -D Trilinos_ENABLE_Shards:BOOL=ON \
  -D Trilinos_ENABLE_ShyLU:BOOL=OFF \
  -D Trilinos_ENABLE_Stokhos:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Teko:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
      -D Teuchos_ENABLE_LONG_LONG_INT:BOOL=ON \
  -D Trilinos_ENABLE_Thyra:BOOL=ON \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Trilinos_ENABLE_TriKota:BOOL=ON \
  -D Trilinos_ENABLE_Zoltan:BOOL=ON \
      -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
  \
  -D Trilinos_ENABLE_Kokkos:BOOL=ON \
  -D Trilinos_ENABLE_KokkosCore:BOOL=ON \
  -D Phalanx_KOKKOS_DEVICE_TYPE:STRING="SERIAL" \
  -D Phalanx_INDEX_SIZE_TYPE:STRING="INT" \
  -D Phalanx_SHOW_DEPRECATED_WARNINGS:BOOL=OFF \
  -D Kokkos_ENABLE_Serial:BOOL=ON \
  -D Kokkos_ENABLE_OpenMP:BOOL=OFF \
  -D Kokkos_ENABLE_Pthread:BOOL=OFF \
  \
  -D Trilinos_ENABLE_STK:BOOL=${HAS_STK} \
      -D Trilinos_ENABLE_STKIO:BOOL=${HAS_STK} \
      -D Trilinos_ENABLE_STKMesh:BOOL=${HAS_STK} \
  \
  -D Trilinos_ENABLE_PyTrilinos:Bool=${HAS_PYTHON} \
  -D PYTHON_EXECUTABLE:PATH=${TACC_PYTHON_BIN}/python3 \
  -D CMAKE_PYTHON_INCLUDE_DIR:PATH="${TACC_PYTHON_INC}" \
  -D CMAKE_PYTHON_LIBRARIES:STRING="${TACC_PYTHON_LIB}" \
  -D PyTrilinos_DOCSTRINGS:BOOL=OFF \
  -D PyTrilinos_ENABLE_Tpetra:BOOL=OFF \
  -D SWIG_EXECUTABLE:FILEPATH=${TACC_SWIG_DIR}/bin/swig \
  -D Trilinos_EXTRA_LD_FLAGS=${PYTHON_LOAD_FLAG} \
  \
  ${srcdir}
#  | tee /admin/build/admin/rpms/stampede2/SPECS/trilinos-${VERSION}-cmake.log 2>&1

export TEMPORARILY_REMOVED="\
  -D VLE_SUPERLU_IS_FIXED_AFTER_12.12:BOOL=ON \
  -D TPL_ENABLE_SuperLU:BOOL=${HAS_SUPERLU} \
      -D SuperLU_INCLUDE_DIRS:PATH="${LMOD_SUPERLUSEQ_DIR}/include" \
      -D SuperLU_LIBRARY_DIRS:PATH="${LMOD_SUPERLUSEQ_LIB}" \
      -D SuperLU_LIBRARY_NAMES:STRING="superlu" \
  \
  "
export albany_extra="\
  -D CMAKE_CXX_FLAGS:STRING=${COPTFLAGS} ${MKLFLAG} -DMPICH_SKIP_MPICXX \
      -DHAVE_AMESOS_SUPERLU5_API -DHAVE_AMESOS2_SUPERLU5_API -DHAVE_IFPACK2_SUPERLU5_API \
  "

# NOX test OFF in an attempt to pass the KNL intel 17 comiler.
# MueLu also OFF because the tutorials don't link

# seems like a bug: https://github.com/trilinos/Trilinos/issues/169

export trilinos_extra_libs="Amesos,Basker,Anasazi,AztecOO,Belos,Epetra,EpetraExt,FEI,Ifpack,Intrepid,ML,MOOCHO,MueLu,NOX,Pamgen,Phalanx,Rhythmos,Sacado,SEACASIoss,SEACAS,SEACASBlot,Shards,ShyLU,Stokhos,Stratimikos,Teko,Teuchos,TriKota,Zoltan; also support enabled for Bool,Hdf5,Netcdf"
