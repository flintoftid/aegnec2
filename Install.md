
# Installation of AEG NEC2

## Configuration

The package uses [CMake][] for configuration. It is strongly recommended to
build the software outside the source tree. Assuming that the source code 
is unpacked into the directory `aegnec2-source` in the current directory 
basic configuration is achieved with the commands

    $ mkdir aegnec2-build
    $ cd aegnec2-build
    $ cmake ../aegnec2-source 

Specific compilers and optimisation flags can be defined using the following syntax:

    $ F77="f77" FFLAGS="-fast -xO4" CC="cc" CFLAGS="-xO2" cmake ../aegnec2-source

This version of [NEC2][] has the ability to use [BLAS][] libraries for the LU
decomposition and back substitution instead of the built in
routines. This may be advantageous if an optimised [BLAS][] library is
available. The implementation uses the [LAPACK][] library from netlib,
however, as a convenience a mini-LAPACK library is shipped with the
aegnec2 source code. By default the configuration script will choose
the solver according to the following rules:

1. If an external LAPACK library is found it will be used.

2. If an external BLAS library is found, but no LAPACK the internal
   mini-LAPACK library will be built and used with the external BLAS.

3. If neither LAPACK or BLAS is found the builtin LU solver is used.

To force the use of a particular solver use

    $ cmake -DSOLVER="BUILTIN" ../aegnec2-source
    $ cmake -DSOLVER="BLAS" ../aegnec2-source
    $ cmake -DSOLVER="LAPACK" ../aegnec2-source

The configuration script will then look for an appropriate library and
fail if it does not locate one.

To force the use of a specific [BLAS][] or [LAPACK][] library, overriding 
the version found by the configuration script use the cmake options
`FIND_BLAS` and `FIND_LAPACK`. Setting these to NO disables the automatic
search for the BLAS and LAPACK libraries which must then be specified
explicitly using `LDFLAGS`, for example:

    $ LDFLAGS=" -lblas " cmake -DSOLVER=BLAS -DFIND_BLAS=NO ../aegnec2 

By default the SOMNEC lossy ground support is built into the `aegnec2`
executable and no external csomnec program is built. This behaviour
can be changed using the following:

    $ cmake -DBUILTIN_SOMNEC=NO../aegnec2-source 

By default installation with be into the `/usr/local` file system. This
can be changed with

    $ cmake -DCMAKE_INSTALL_PREFIX=/somepath ../aegnec2-source 

where `/somepath` will be prefixed to the installation paths.


## Compilation

Type

    make

to build the executables. The package contains a regression test suite 
which can be run using

    make test

This assumes enough memory is available to run a jobs with a 64MB core
size (2000 segments). If you do not have this much memory try deleting
the file tests/test2000.nec before running the test suite. This
will limit the test suite to requiring around 24MB of RAM.

If any tests are failed check the appropriate log file in the
test suite directory to see if the problem is serious. The validation
assumes that the currents and fields calculated by `aegnec2` are identical
to those on the the reference system. This requires, at least, that 
your system is IEEE 754 floating point compliant. Slight variation in 
the last couple of decimal places in currents and fields may occur under
certain circumstances on some systems.

## Installation

The `aegnec2` and `somnec` executables and the man pages can be installed
with

    make install

By default this will place them in `/usr/local/bin`. See above for 
instructions on changing the installation prefix. Alternatively
a binary package can be built using

    make package

and then installed separately. 

### Installation of Binary Packages

Two types of binary package are available: a compessed tar file and 
a self-extracting Bourne shell installer.


## Platform specific notes

### Linux x86 & GNU Compilers

aegnec2 can be successfully built using [GCC][] 3.4 or earlier and
gcc/gfortran 4.3 or above. Earlier versions of gfortran may have
some problems.

### Linux x86 & Portland Group Inc. Compilers

Using the Portland Group compilers [PGI][] below version 3.0 on Linux the
`aegnec2` executable fails to link (due to contention between C and F77
over who controls the entry point). The executable can be link
manually by changing into the cnec directory and doing something like:

    pgcc  -O2 -tp p6 -Munroll -Mdalign -L../lib -o aegnec2  nec.o cnec2d.o \
          ../nec2d/libnec2d.a ../somnec/libsomnec.a ../somnec/libsomtest.a \
          ../mlapack/libmlapack.a -lgnu -L/scratch/idf1/lib -lblas \
          -L/scratch/idf1/pgi/linux86/lib -lpgftnrtl

where the PGI compiler is assumed to be installed in `/scratch/pgi`. The
key point is to use the C compiler to link the executable and link
`libpgftnrtl.a` but not `pgfmain.o` from the PGI `lib` directory.  For
version 3.0 and above of the PGI compilers this is solved
transparently by the new `-Mnomain` flag.

### Linx x86, Intel Compilers and Math Kernel Library

To use the buitl-in mini LAPACK:

    $ CC=icc FC=ifort LDFLAGS="-L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32" cmake -DSOLVER=BLAS -DFIND_BLAS=NO ../aegnec2 

To use the Intel LAPACK:

    $ CC=icc FC=ifort LDFLAGS="-L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32" cmake -DSOLVER=LAPACK -DFIND_LAPACK=NO ../aegnec2

### Microsoft Windows & Cygwin
 
The Cygnus [Cygwin][] GNU development environment can be used to
compile aegnec2 with either the builtin solver or using the reference
[BLAS][] libraries from netlib. Instructions on installing and using 
these packages can be found in their own documentation.

To compile executables for the [Cygwin][] environment the standard
Unix procedure is adequate.


[NEC2]:                      https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code
[Mercurial]:                 http://mercurial.selenic.com
[CMake]:                     http://www.cmake.org
[C89]:                       https://en.wikipedia.org/wiki/ANSI_C#C89
[Fortran77]:                 https://en.wikipedia.org/wiki/Fortran
[LAPACK]:                    http://www.netlib.org/lapack
[BLAS]:                      http://www.netlib.org/blas
[GCC]:                       https://gcc.gnu.org
[PGI]:                       http://www.pgroup.com
[Cygwin]:                    http://www.cygwin.com
[MINGW]:                     http://www.mingw.org
[ATLAS]:                     http://math-atlas.sourceforge.net
[MKL]:                       http://software.intel.com/en-us/intel-mkl
