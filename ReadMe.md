![](https://github.com/flintoftid/aegnec2/blob/master/doc/aegnec2.jpg)

# AEG NEC2: An open source MoM solver for electromagnetic simulations

The Applied Electromagnetics Group ([AEG][]) boundary element method ([BEM][]) 
method-of-moments ([MoM][]) solver, `aegnec2`, is an [Open Source][] derivative of 
the Numerical Electromagnetics Code Version 2 ([NEC2][]). It was developed in the in 
the [Department of Electronic Engineering][] at the [University of York][] for research in 
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

## Documentation

Installation instructions are contained in the file [Install.md][] in the 
source distribution.

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

Please report bugs using the gitgub issue tracker at
<https://github.com/flintoftid/aegnec2/issues> or by email to <ian.flintoft@googlemail.com>.

For general guidance on how to write a good bug report see, for example:

* <http://www.chiark.greenend.org.uk/~sgtatham/bugs.html>
* <http://noverse.com/blog/2012/06/how-to-write-a-good-bug-report>
* <http://www.softwaretestinghelp.com/how-to-write-good-bug-report>

Some of the tips in <http://www.catb.org/esr/faqs/smart-questions.html> are also relevant to reporting bugs.

## How to contribute

We welcome any contributions to the development of the code, including:

* Fixing bugs.

* Interesting examples that can be used for test-cases.

* Improving the user documentation.

Please contact [Dr Ian Flintoft][], <ian.flintoft@googlemail.com>, if you are 
interested in helping with these or any other aspect of development.

## Licence

The code is licensed under the GNU Public Licence, version 3 ([GPL3][]). 

The mlapack directory contains a few routines from the 
[LAPACK 3.2.1 mathematics libraries](http://www.netlib.org/lapack) which are
built into a mini LAPACK library under certain build options. 
These files are distributed under their own GPL compatible license.

See the file [Licence.md][] in the souce code distribution for details.

## Developer 

[Dr Ian Flintoft][], <ian.flintoft@googlemail.com>

## Contacts

[Dr Ian Flintoft][] : <ian.flintoft@googlemail.com>

[Dr John Dawson][] : <john.dawson@york.ac.uk>

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

The source code in the lapack directory is from [netlib][]. See individual
source files for author information.

For information on the various source files see [Licence.md][] in the source or binary package. 
For a history of NEC codes see ([Adler1993][]).

## Publications using AEG NEC2

[Capstick2009]: http://dx.doi.org/10.1109/TIM.2008.2005263

([Capstick2009]) M. H. Capstick, J. Jekkonen, A. C. Marvin, I. D. Flintoft and L. Dawson, 
“A novel indirect method to determine the radiation impedance of a handheld antenna structure”, 
IEEE Transactions on Instrumentation and Measurement, vol. 58, no. 3, pp. 578-585, 2009.

[Jekkonen2005b]: http://www.ursi.org/publications.php

([Jekkonen2005b]) J. Jekkonen, I. D. Flintoft, M. H. Capstick and A. C. Marvin, 
“Application of indirect measurements techniques to mobile communication antenna systems”,
XXVIIIth Triennial General Assembly of the International Union of Radio Science, New Delhi, 
India, paper no. B05P.11(072), 23-29 Oct. 2005.

(Jekkonen2005a) J. Jekkonen, I. D. Flintoft, M. H. Capstick and A. C. Marvin, “A novel 
indirect method to determine the radiation impedance of an unknown antenna structure”, 
16th International Zurich Symposium and Technical Exhibition on Electromagnetic Compatibility, 
pp. 179-182, Zurich, Switzerland, 14-18 Feb. 2005.

[Papatsoris2004]: http://dx.doi.org/10.1049/ip-smt:20040566

([Papatsoris2004])	A. D. Papatsoris, I. D. Flintoft, D. W. Welsh and A. C. Marvin, 
“Modelling the cumulative emission field of unstructured telecommunication 
transmission networks”, IEE Proceedings on Science, Measurement and Technology, 
vol. 151, no. 4, pp. 244-252, 2004.

(Flintoft2003) I. D. Flintoft, A. D. Papatsoris, D. W. Welsh, A. C. Marvin, “Radiated 
emissions from unstructured networks: Potential impact on maritime and aeronautical radio 
services”, 15th International Zurich Symposium and Technical Exhibition on Electromagnetic 
Compatibility, Zurich, Switzerland, pp. 93-98, 18-20 Feb. 2003.

[Papatsoris2002]: http://dx.doi.org/10.1049/el:20021090

([Papatsoris2002]) A. D. Papatsoris and I. D. Flintoft, “Model for the emission 
electric field of distributed unstructured telecommunication transmission networks”, 
IEE Electronics Letters, vol. 38, no. 24, pp.1610-1611, 2002.

[Jekkonen2002]: http://www.ursi.org/publications.php

([Jekkonen2002]) J. Jekkonen, A. C. Marvin and I. D. Flintoft, “An indirect method for 
measuring the radiation impedance of an unknown antenna structure”, XXVIIth Triennial General 
Assembly of the International Union of Radio Science, Maastricht, The Netherlands, 
paper no. AE.O.4, 17-24 Aug. 2002.

(Welsh2001) D. W. Welsh, A. D. Papatsoris, I. D. Flintoft and A. C. Marvin, “Investigation 
of likely increases in established radio noise floor due to widespread deployment of PLT, 
ADSL and VDSL broadband access technologies”, 14th International Zurich Symposium and 
Technical Exhibition on Electromagnetic Compatibility, Zurich, Switzerland, pp. 595-600, 
20-22 Feb 2001.

[Papatsoris2000]: http://dx.doi.org/10.1049/el:20000841

([Papatsoris2000]) A. D. Papatsoris and I. D. Flintoft, “Modelling the cumulative 
ground wave electric fields from the widespread deployment of xDSL data distribution 
systems”, IEE Electronics Letters, vol. 36, no. 13, pp. 1171-1172, 2000.

(Flintoft1998) I. D. Flintoft, S. J. Porter and A. C. Marvin, “Interaction of wired IT 
networks and mobile telecommunication systems”, 3rd European Symposium on Electromagnetic
Compatibility (EMC’98 Roma), Rome, Italy, pp. 832-836, 14-18 Sep. 1998.

## Related links

* Other versions of [NEC2][] can be found in the [Unofficial NEC Archives][].

* The original [NEC2 Report][], G J Burke and A J Poggio, Numerical Electromagnetics Code (NEC), 
  Technical Report UCID-18834, Lawrence Livermore National Laboratories, 1981.

## References

[Adler1993]: http://www.aces-society.org/search.php?vol=8&no=3&type=3

([Adler1993]) D. Adler, Software Exchange Committee Report, ACES Newsletter, vol. 8, no. 3, pp. 8-10, November 1993.


[NEC2]: https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code
[NEC2 Report]: http://www.dtic.mil/dtic/tr/fulltext/u2/a956129.pdf
[NEC2 Manuals]: http://www.nec2.org
[Unofficial NEC Archives]: http://nec-archives.pa3kj.com
[GPL3]: http://www.gnu.org/copyleft/gpl.html
[Open Source]: http://opensource.org
[CMake]: http://www.cmake.org
[C89]: https://en.wikipedia.org/wiki/ANSI_C#C89
[Fortran77]: https://en.wikipedia.org/wiki/Fortran
[netlib]: http://www.netlib.org
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
[Department of Electronic Engineering]: https://www.york.ac.uk/electronic-engineering
[AEG]: https://www.york.ac.uk/electronic-engineering/research/communication-technologies/applied-electromagnetics-devices
[Dr Ian Flintoft]: https://flintoftid.github.io
[Dr John Dawson]: https://www.york.ac.uk/electronic-engineering/staff/john_dawson
[Dr Stuart Porter]: https://www.york.ac.uk/electronic-engineering/staff/stuart_porter
[EMC]: https://www.york.ac.uk/electronic-engineering/research/communication-technologies/applied-electromagnetics-devices/emc-shielding/
[CEM]: https://www.york.ac.uk/electronic-engineering/research/communication-technologies/applied-electromagnetics-devices/numerical-modelling-optimisation/
[Licence.md]: https://github.com/flintoftid/aegnec2/blob/master/Licence.md
[Install.md]: https://github.com/flintoftid/aegnec2/blob/master/Install.md
[Bugs.md]: https://github.com/flintoftid/aegnec2/blob/master/doc/Bugs.md
[ToDo.md]: https://github.com/flintoftid/aegnec2/blob/master/doc/ToDo.md

