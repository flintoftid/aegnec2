![](https://bitbucket.org/uoyaeg/aegnec2/wiki/image.jpg )

aegnec2: An open source MoM solver for electromagnetic simulations
==================================================================

Introduction
------------

*aegnec2* is a derivative of the Numerical Electromagnetics Code Version 2
([NEC2][]) with a number of enhancements including:

* Dynamically allocated memory using command line switches.

* Ability to use [BLAS][] based LU solvers from [LAPACK][] for increased
  performance if available.

* Fully double precision including all embedded constants.

* Improved user interface.

* Regression test suite.

* The `somnec` lossy ground model can be built directly into the nec
  binary making lossy ground calculations more straightforward and
  more flexible.

 * Bug fixes and portability improvements.

For authors and copyright information on the various source files see
the files [Authors.txt][] and [Licence.txt][] in the source or binary 
package. For a history of NEC codes see [ACES1993][].

Requirements and portability
----------------------------

Compilation requires a [Fortran 77][] compiler with the following
extensions:

* Preprocessing of `.F` files either natively or transparently via `cpp`.

* `INCLUDE` statement to include files.

* `REAL*8` and `COMPLEX*16` data types.

* The following intrinsic functions: `DREAL`, `DIMAG`, `DCMPLX`, `DCONJG`,
  `CDABS`, `CDSQRT`, `CDEXP` and `DFLOAT`.

* List directed I/O from internal files.

* Portable equivalencing of `INTEGER` and `REAL*8`. Explicitly
  this must work:

    INTEGER I1 , I2
    REAL*8  R1

    COMMON /FRED/ I1 ,  I2

    EQUIVALENCE (I1,R1)

  such that R1 occupies the same storage as two INTEGER's. This
  restriction is likely to removed in a later version.

Timing information requires that the Fortran compiler supports one of
the following intrinsics: `ETIME`, `CPU_TIME` or `SECOND`. Otherwise all times 
will be reported as zero.

An ANSI [C89][] compiler is also required for the front end. The C to
Fortran interface is currently very basic. It requires that the C data
types and Fortran data types map like:

C type   | Fortran type
:--------|:-------------
`int`    | `INTEGER`
`double` | `REAL*8`

If also requires that the length of character strings is passed by value 
as an `unsigned long` at the end of the argument list. Most of this can 
easily be generalised to work on other systems if necessary. Please provide
patches to the file `nec/cnec2d.c` only if possible to generalise the interface. 

A level 3 [BLAS][] library is required to use the [LAPACK][] based LU
solver. The reference implementation is available from the netlib
archives [BLAS][]. Many OS and compiler vendors supply highly optimised
versions of the BLAS libraries. Optimised BLAS/LAPACK libraries that 
have been tested with aegnec2 include:

 * The ASCI Red Pentium Pro Optimised BLAS Libraries developed by Greg
   Henry for GNU/Linux on Pentium Pro's. PII's and PIII's [ASCI Red][].

 * The portable [ATLAS][] library which provides a code generator for
   creating optimised BLAS (and some LAPACK) routines on a variety of
   platforms.
 
 * The Intel Math Kernel Library ([MKL][]) for GNU/Linux.

For development of the code [CMake][] version 2.6.4 or higher is required.

Installation from source on Linux/Unix/Cygwin
---------------------------------------------

### Configuration

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

### Compilation

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

### Installation

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

### Platform Specific Notes

#### Linux x86 & GNU Compilers

aegnec2 can be successfully built using [GCC][] 3.4 or earlier and
gcc/gfortran 4.3 or above. Earlier versions of gfortran may have
some problems.

#### Linux x86 & Portland Group Inc. Compilers

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

#### Linx x86, Intel Compilers and Math Kernel Library

To use the buitl-in mini LAPACK:

    $ CC=icc FC=ifort LDFLAGS="-L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32" cmake -DSOLVER=BLAS -DFIND_BLAS=NO ../aegnec2 

To use the Intel LAPACK:

    $ CC=icc FC=ifort LDFLAGS="-L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32" cmake -DSOLVER=LAPACK -DFIND_LAPACK=NO ../aegnec2

#### Microsoft Windows & Cygwin
 
The Cygnus [Cygwin][] GNU development environment can be used to
compile aegnec2 with either the builtin solver or using the reference
[BLAS][] libraries from netlib. Instructions on installing and using 
these packages can be found in their own documentation.

To compile executables for the [Cygwin][] environment the standard
Unix procedure is adequate.


Installation From Source on Microsoft Windows
---------------------------------------------

Not supported yet!

Troubleshooting
---------------

### aegnec2 segmentation faults almost instantly when running my model. What is wrong?

You may have specified fewer segments or patches than your job actually requires. Try 
using the `-s` and/or `-p` options to increase that array sizes and see if works then.


Using aegnec2
-------------

Basic usage information is available using the `-h` or `--help` flags:

    aegnec2 -h

or from the man pages. The input file follows the same syntax as
"standard" NEC2 except that full line and end of line comments are
allowed using a hash (#). This is useful for temporarily commenting
out a few cards. The input file name must end with a `.nec`
suffix. This may be relaxed in a later version. The suffix can be
omitted when running `aegnec2`. On-line [NEC2 manuals][] are available
and there is a [NEC Mailing List][]. Other versions of the code can
be found in the [Unofficial NEC Archives][].

If `aegcnec2` is run with no options the default behaviour is similar to a
300 unknown version of "standard" NEC2. If you job has more than 300
segments you will need to specify the number with the -s option, for
example, if your job has 2000 segments then

    aegnec2 -s 2000 job.nec

will run NEC2 with sufficient storage for 2000 segments. You can
specify a FEW more segments than you actually need, however, it is
anti-social to run jobs which hog lots of memory for no purpose.

By default the core size of your job will be matched to exactly fit
the number of unknowns you specify. It is recommended that you run
`aegnec2` in this way.  However for very large jobs running on machines
with limited memory or many users it may be beneficial to limit the
core size. This can be done with the `-c` option. For example to run a
2000 segment job, which normally would use 64Mb of core storage, with
only 48Mb of core use

    aegnec2 -s 2000 -c 48m job.nec

This will use an out of core solver which is generally slower than an
in core solution but may be faster if the job would otherwise cause
your machine to use significant swap space. The suffices b, k and m
are used to denote bytes, kilobytes (1024 bytes) and megabytes (1024
kilobytes) respectively.

By default all the results are written into a file call job.res
(assuming the input file is called `job.nec`). If you use the `-M`
(`--old-multi-file`) option, e.g.

    aegnec2 -M job.nec

then any near field and far field data will also be written to
separate files called

File      | Contents
:---------|:-----------------------------------
`job.nre` | Near electric field data (NE card)
`job.nrh` | Near magnetic field data (NH card)
`job.rdp` | Radiation patterns (RP card)

for easier plotting. Impedance data will also be written to job.imp.
If multiple NE, NH or RP cards are used ALL the data is written to the
same file in the order in which it is calculated. If you would prefer
each NE, NH and RP card to put it's data in a separate file then use
the -m option (`--multi-file`):

    aegnec2 -m job.nec

In this case a separate file is created each time a card is used,
e.g. if four NE cards are used (or one NE with four frequency points):

    job-0001.nre
    job-0002.nre
    job-0003.nre
    job-0004.nre

will be created.

References
----------

[ACES1993]:                  http://www.aces-society.org/newsletter.php          
                             "ACES Newsletter, vol. 8, no. 3, pp. 8-10, November 1993"
[NEC2]:                      http://www.dtic.mil/dtic/tr/fulltext/u2/a956129.pdf 
                             "G J Burke and A J Poggio, Numerical Electromagnetics Code (NEC), 
                             Technical Report UCID-18834, Lawrence Livermore National Laboratories, 1981"
[NEC2 Manuals]:              http://www.nec2.org                                 
                             "On-line NEC2 Manuals"
[NEC Mailing List]:          http://www.robomod.net/mailman/listinfo/nec-list    
                             "NEC mailing list archives"
[Unofficial NEC Archives]:   http://nec-archives.pa3kj.com/                      
                             "Unofficial Numerical Electromagnetic Code (NEC) Archives"

[Open Source]:               http://opensource.org
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

[University of York]:        http://www.york.ac.uk
[Department of Electronics]: http://www.elec.york.ac.uk
[AEG]:                       http://www.elec.york.ac.uk/research/physLayer/appliedEM.html
[Dr Ian Flintoft]:           http://www.elec.york.ac.uk/staff/idf1.html
[Dr John Dawson]:            http://www.elec.york.ac.uk/staff/jfd1.html
[Dr Stuart Porter]:          http://www.elec.york.ac.uk/staff/sjp1.html

[Install.txt]:               https://bitbucket.org/uoyaeg/aegnec2/raw/tip/Install.txt
[Bugs.txt]:                  https://bitbucket.org/uoyaeg/aegnec2/raw/tip/doc/Bugs.txt
[ToDo.txt]:                  https://bitbucket.org/uoyaeg/aegnec2/raw/tip/doc/ToDo.txt
