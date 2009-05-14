#!/bin/sh

ARCH_FLAGS=" -tp p6 "
OPT_FLAGS=" -O2 -Munroll -Mdalign "

CC=pgcc CFLAGS="$OPT_FLAGS $ARCH_FLAGS" FC=pgf77 FFLAGS="$OPT_FLAGS $ARCH_FLAGS" cmake -DCPACK_SYSTEM_NAME="linux-i686" $*
make 
make package
make package_source

exit $?

