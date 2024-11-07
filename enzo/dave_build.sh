#!/bin/bash

cd ${STOCKYARD}/enzo

config_common="-DCHARM_ROOT="${TACC_CHARMPP_DIR}" \
    -Duse_projections=OFF \
    -DUSE_YT_BASED_TESTS=OFF \
    -DBoost_INCLUDE_DIR=$HOME/include \
    -DBoost_LIB_DIR=${TACC_BOOST_LIB} \
    -DUSE_GRACKLE=OFF \
    -DUSE_DOUBLE_PREC=OFF \
    -DBUILD_TESTING=ON \
    -DNEW_FIELD=ON \
    -Duse_jemalloc=OFF \
    -Dbalance=OFF  \
    -Ddarshan=OFF \
    -Dmemory=OFF \
    -Dsmp=OFF \
 -DGRACKLE_INPUT_DATA_DIR=$HOME/data/grackle/grackle_data_files/input"

# -DGrackle_ROOT=$HOME \
# -DPNG_PNG_INCLUDE_DIR=$HOME/include \
# -DPNG_LIBRARY=$HOME/lib/libpng.so \
#

rm -rf dave_build
mkdir dave_build
cd dave_build

cmake $config_common -DEnzo-E_CONFIG=linux_nvc ../enzo-git
make

