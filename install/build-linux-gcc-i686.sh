#!/bin/sh

ARCH_FLAGS="-mfpmath=sse -malign-double -march=pentium4"
OPT_FLAGS="-O2"

CC=gcc CFLAGS="$OPT_FLAGS $ARCH_FLAGS" FC=gfortran FFLAGS="$OPT_FLAGS $ARCH_FLAGS" cmake -DCPACK_SYSTEM_NAME="linux-i686" $*
make 
make package
make package_source

exit $?

