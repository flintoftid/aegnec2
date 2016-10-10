![](https://bitbucket.org/uoyaeg/aegnec2/wiki/aegnec2.jpg )

# AEG NEC2: An open source MoM solver for electromagnetic simulations

The Applied Electromagnetics Group ([AEG][]) boundary element method ([BEM][]) 
method-of-moments ([MoM][]) solver, `aegnec2`, is an [Open Source][] derivative of 
the Numerical Electromagnetics Code Version 2 ([NEC2][]). It was developed in the in 
the [Department of Electronics][] at the [University of York][] for research in 
electromagnetic compatibility ([EMC][]) and computational electromagnetics ([CEM][]).


## Code Features

`aegnec2` has a number of enhancements compared to the original version of ([NEC2][]) including:

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


## Requirements

### CMake (mandatory)

In order to compile and install the code from source using the build system provided a 
recent version of [CMake][] is required.

### Fortran77 compiler (mandatory)

The core code is written in [Fortran77][] and requires a compiler with the 
following [Fortran77][] extensions:

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

### ANSI C compiler (mandatory)
 
An ANSI [C89][] compiler is required for the front end. The C to
Fortran interface uses CMakes built-in support. It requires that the C data
types and Fortran data types map like:

C type   | Fortran type
:--------|:-------------
`int`    | `INTEGER`
`double` | `REAL*8`

It also requires that the length of character strings is passed by value 
as an `unsigned long` at the end of the argument list.

### BLAS libraries (optional)

A level 3 [BLAS][] library is required to use the [LAPACK][] based LU
solver. The reference implementation is available from the netlib
archives [BLAS][]. Many OS and compiler vendors supply highly optimised
versions of the BLAS libraries. Optimised BLAS/LAPACK libraries that 
have been tested with aegnec2 include:

* The portable [ATLAS][] library which provides a code generator for
  creating optimised BLAS (and some LAPACK) routines on a variety of
  platforms.
 
* The Intel Math Kernel Library ([MKL][]) for GNU/Linux.


## Installation

Installation instructions are contained in the file [Install.md][] in the 
source distribution.


## Documentation

Basic usage information for the `aegnec2` and `aegsomnec2` programs is 
available from the executables using the `-h` or `--help` flags, e.g.

    $ aegnec2 -h

or from the man pages

    $ man aegnec2
    
after installation . The input file follows the same syntax as
"standard" NEC2 except that full line and end of line comments are
allowed using a hash (`#`). This is useful for temporarily commenting
out a few cards. The input file name must end with a `.nec`
suffix. This may be relaxed in a later version. The suffix can be
omitted when running `aegnec2`. On-line [NEC2 manuals][] are available.

If `aegcnec2` is run with no options the default behaviour is similar to a
300 unknown version of "standard" NEC2. If you job has more than 300
segments you will need to specify the number with the -s option, for
example, if your job has 2000 segments then

    $ aegnec2 -s 2000 job.nec

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

    $ aegnec2 -s 2000 -c 48m job.nec

This will use an out of core solver which is generally slower than an
in core solution but may be faster if the job would otherwise cause
your machine to use significant swap space. The suffices b, k and m
are used to denote bytes, kilobytes (1024 bytes) and megabytes (1024
kilobytes) respectively.

By default all the results are written into a file call job.res
(assuming the input file is called `job.nec`). If you use the `-M`
(`--old-multi-file`) option, e.g.

    $ aegnec2 -M job.nec

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

    $ aegnec2 -m job.nec

In this case a separate file is created each time a card is used,
e.g. if four NE cards are used (or one NE with four frequency points):

    job-0001.nre
    job-0002.nre
    job-0003.nre
    job-0004.nre

will be created.


## Bugs and support

Known significant bugs are listed in the file doc/[Bugs.md][] in the source code. 

Please report bugs using the bitbucket issue tracker at
<https://bitbucket.org/uoyaeg/aegnec2/issues> or by email to <ian.flintoft@googlemail.com>.

For general guidance on how to write a good bug report see, for example:

* <http://www.chiark.greenend.org.uk/~sgtatham/bugs.html>
* <http://noverse.com/blog/2012/06/how-to-write-a-good-bug-report>
* <http://www.softwaretestinghelp.com/how-to-write-good-bug-report>

Some of the tips in <http://www.catb.org/esr/faqs/smart-questions.html> are also relevant to reporting bugs.


## Licence

The code is licensed under the GNU Public Licence, version 3 ([GPL3][]). 

The mlapack directory contains a few routines from the 
[LAPACK 3.2.1 mathematics libraries](http://www.netlib.org/lapack) which are
built into a mini LAPACK library under certain build options. 
These files are distributed under their own GPL compatible license.

See the file [Licence.md][] in the souce code distribution for details.


## Developer 

[Dr Ian Flintoft][], <ian.flintoft@googlemail.com>


## Credits

Many people have contributed to NEC2 including:

R. W. Adams
J. N. Brittingham
G. J. Burke
F. J. Deadrick
K. K. Hazard
D. L. Knepp
D. L. Lager
R. J. Lytle
E. K. Miller
J. B. Morton
G. M. Pjerrou
A. J. Poggio
E. S. Selden

This version of NEC2, called aegnec2, is maintained by [Dr Ian Flintoft][], <ian.flintoft@googlemail.com>.

The source code in the lapack directory is from netlib. See individual
source files for author information.

For information on the various source files see [Licence.md][] in the source or binary package. 
For a history of NEC codes see ([Adler1993][]).


## Related links

* The [NEC Mailing List][]. 

* Other versions of [NEC2][] can be found in the [Unofficial NEC Archives][].

* The original [NEC2 Report][], G J Burke and A J Poggio, Numerical Electromagnetics Code (NEC), 
  Technical Report UCID-18834, Lawrence Livermore National Laboratories, 1981.


## References

[Adler1993]: http://www.aces-society.org/search.php?vol=8&no=3&type=3

([Adler1993]) D. Adler, Software Exchange Committee Report, ACES Newsletter, vol. 8, no. 3, pp. 8-10, November 1993.


[NEC2]: https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code
[NEC2 Report]: http://www.dtic.mil/dtic/tr/fulltext/u2/a956129.pdf
[NEC2 Manuals]: http://www.nec2.org
[NEC Mailing List]: http://www.robomod.net/mailman/listinfo/nec-list
[Unofficial NEC Archives]: http://nec-archives.pa3kj.com

[GPL3]: http://www.gnu.org/copyleft/gpl.html
[Open Source]: http://opensource.org
[Mercurial]: http://mercurial.selenic.com
[CMake]: http://www.cmake.org
[C89]: https://en.wikipedia.org/wiki/ANSI_C#C89
[Fortran77]: https://en.wikipedia.org/wiki/Fortran
[LAPACK]: http://www.netlib.org/lapack
[BLAS]: http://www.netlib.org/blas
[GCC]: https://gcc.gnu.org
[PGI]: http://www.pgroup.com
[Cygwin]: http://www.cygwin.com
[MINGW]: http://www.mingw.org
[ATLAS]: http://math-atlas.sourceforge.net
[MKL]: http://software.intel.com/en-us/intel-mkl
[MoM]: https://en.wikipedia.org/wiki/Computational_electromagnetics
[BEM]: https://en.wikipedia.org/wiki/Boundary_element_method
[University of York]: http://www.york.ac.uk
[Department of Electronics]: http://www.elec.york.ac.uk
[AEG]: http://www.elec.york.ac.uk/research/physLayer/appliedEM.html
[Dr Ian Flintoft]: https://idflintoft.bitbucket.io
[Dr John Dawson]: http://www.elec.york.ac.uk/staff/jfd1.html
[Dr Stuart Porter]: http://www.elec.york.ac.uk/staff/sjp1.html
[EMC]: http://www.elec.york.ac.uk/research/physLayer/appliedEM/emc.html 
[CEM]: http://www.elec.york.ac.uk/research/physLayer/appliedEM/numerical.html

[Licence.md]: https://bitbucket.org/uoyaeg/aegnec2/raw/tip/Licence.md
[Install.md]: https://bitbucket.org/uoyaeg/aegnec2/raw/tip/Install.md
[Bugs.md]: https://bitbucket.org/uoyaeg/aegnec2/raw/tip/doc/Bugs.md
[ToDo.md]: https://bitbucket.org/uoyaeg/aegnec2/raw/tip/doc/ToDo.md

