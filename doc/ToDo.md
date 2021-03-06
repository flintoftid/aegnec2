
# AEG NEC2 To Do List


## Version 0.8.1

### Release objectives for 0.8.1

1. Initial Windows support via MINGW and NSIS.

2. Performance assessment with MINGW.

3. Documentation for building with MINGW.

4. Tested (not via test suite) on WinXP/Intel.

### Help screen won't fit into a DOS box in Win95/98/NT. 

Need an abbreviated help or some way to work around the problem. 

    cnec2 -h|more
    
doesn't seem to work. Maybe short help with `-h` (note to use `-V` as well 
for more help) and long help with `-V -h` and a `"Hit return to continue"` 
set up. How many lines per screen?

### Add more support for building under Windows

Including using the Intel Maths Kernel Libraries. How do you do this? 
Make a windows source distribution using some additional targets in the 
make file. Could for example combine all the Fortran code into a single 
file and pre-process for known windows support. May need to leave HAVE_BLAS
as a definable. 

### Investigate how to get optimised BLAS in Win32. 

1. Cygwin/mingw: Either need source (asm) for an optimised BLAS or a way
   of using the Intel Performance Library DLL with mingw. Cygwin is slow
   anyway.

2. Native MS compilers: should work with Intel Performance Library but need
   a "make" system for compilation. Is there a portable way to do this for
   all Win32 compilers/IDEs etc (dream on)?

3. ATLAS.


## Version 0.8.2

### Release objectives for 0.8.2

1. Port testsuite to work on Win32 and Unix.

### Fix broken tests

Issues are negative zero's: 0.00 versus -0.00. Could preprocess results 
file and change any negative zero to positive. Could be different number 
of digits after decimal point. Need some serious regexp!

### Convert testsuite drivers to Python

### More robust way to extract timing data

### More robust way to perform comparision of results file

The validation of `testngf1.nec` and `testngf2.nec` is known to fail
on alpha platforms, both Tru64 Unix and Linux. This is due to a
precision problems in regenerating the geometry data from the binary
NGF file causing `tan(0,x)` to be 180 instead of 0 due to noise in small
values of `x`. `atgn2.f` wraps the standard `DATAN2(Y,X)` function such that
it returns `0` if `X==Y==0`. If `Y==0` but `X` is a very small number then
rounding in the NGF code can flip its sign causing `ATGN2` to flip from
`0` to `180`. These values are effectively degenerate so the failure is
not significant. It is however annoying and will be fixed eventually.
Can this be solved using a tolerance? `"DX" <= 1D-15` goes to zero? Can 
the calculation in GFIL be improved for small vales of `X, Y ,Z`? Could 
this be done in the atgn2 routine or is this too low level? For 
example (note changed XX<-->YY for sanity):
 
    REAL*8 FUNCTION ATGN2 ( YY , XX )
    IMPLICIT NONE
    REAL*8 XX , YY
    IF( XX.LT.1.0D-15 ) THEN
       IF( YY.LT.1.0D-15 ) THEN
          ATGN2 = 0.0D0
       ELSE
          ATGN2 = DATAN2(YY,XX)	
       ENDIF
    ELSE
       ATGN2 = DATAN2(YY,XX)	
    ENDIF
    RETURN
    END

Does this incur a significant performance penalty? Do I need test on
YY in outer ELSE statement?


## Future versions

### Run sphere model again on new machine

### Add code to wire.f to catch adding wires beyond end of array sizes (LD)

Presumably there will be other places where this may occur - candidates are:

    emc5 $ grep "^ *X(.*) *=" *.f *.F
    arc.f:         X(I) = XS1
    conect.f:            X(IX) = XA
    datagn.f:         X(I) = X(I)*XW1
    datagn.f:         X(I) = X(I)*XW1
    datagn.f:         X(I) = (X(I)+X2(I))*0.5D0
    gfil.f:         X(I) = XI - ALP(I)*DX
    gfil.f:         X(I) = X(I)*WLAM
    helix.f:         X(I) = A1*DCOS(2.0D0*PI*Z(I)/SSS)
    helix.f:         X(I) = (A1+(A2-A1)*Z(I)/DABS(HL))*DCOS(2.0D0*PI*Z(I)/SSS)
    helix.f:         X(I) = Y(I)
    move.f:            X(K) = XI*XX + YI*XY + ZI*XZ + XXS
    move.f:            X(KR) = XI*XX + YI*XY + ZI*XZ + XXS
    patch.f:      X(MI) = X1
    patch.f:      X(MI) = X1 + 0.5D0*(S1X+S2X)
    patch.f:      X(MI) = (X1+XXX2+X3)/3.0D0
    patch.f:      X(MI) = (XA*(X1+XXX2+X3)+XST*(X1+X3+X4))*SALPN
    patch.f:            X(MI) = XN2 + XST*S1X
    reflc.f:         X(I) = XK*CS - YK*SS
    reflc.f:         X(J) = XK*CS - YK*SS
    subph.f:         X(NYP) = X(IX)
    subph.f:         X(MIA) = XXS + XT*S1X + YT*S2X
    wire.f:         X(I) = XS1
    nec2d.F:         X(I) = XTEMP(I)*FR
    nec2d.F:         X(J) = XTEMP(J)*FR
 
### Add documentation for PL card

    PL IPLP1 IPLP2 IPLP3 IPLP4

    IPLP1 2 Near pattern printed to channel 8 (CHRPAT)
    IPLP2 1 Print fields as complex numbers i.e. real and imag
          2 Print fields as mag and phase
    IPLP3 1 Print EX only
          2 Print EY only
          3 Print EZ only
          4 Print EX, EY and EZ
          5 Print ETOT (IPLP2=2 only)
    IPLP4 1 Print X pos of obs point
          2 Print Y pos of obs point
          3 Print Z pos of obs point

`datagn.f`:

    IPLP1=1 
	  
PL also seems to have effects in `datagn.f`, `nec2.F`, `rdpat.f` and
possibly else where - check these out and document.

    PL 1 1 1 1

Currents on segment in real and imaginary form

### See if possible to remove equivalenced INTEGER and REAL*8 types

If so remove configure test as well. Do in both C interface and `fdrive.f`.
What is impact on memory usage of doing this (as a function of core
size)? Want a graph of percentage memory increase versus number of
segments. All the equivalenced arrays are at most O(N) so effect should 
rapidly decline with N.

### Consider removing preprocessed Fortran code

Can do this using autoconf to process `timer.f`.in for example or substituting 
source files, eg `timer-etime.f`, `timer-cpu_time.f` etc using autoconf. For
timing stuff later option is probably better. For blas support this 
 option would really have to be used (`factrs-nec.f`, `factrs-blas.f`, 
`solves-nec.f`, `solves-blas.f`) since a number of lines of code are 
involved.

Could `INCLUDE` be used to remove dependence of a preprocessor?
E.g FLUSH case: May also be able to do generated INCLUDE `flush.inc` 
file which is blank if not supported.

### See if constants in ZINT can be recalculated to greater precision using MuPAD

Yes - fairly straightforward (A+S). Need more terms in
series? (#ifdef?). Details are in my LINK Project Note Book (don't
know when I did this)! See scanned file `kelvin.pdf`.

### Remake headers with gxchk

### PATOFF needs to be greater than or equal to NUM_SEG_MAX in nec.h

Currently both are set to 10000 (standard NEC2 offset) for
safety. Need to get this limit from the same place in both C and
Fortran. What is the 10000.0D0 in subph for?

### File suffices should be set in just one place

Also all buffers should automatically scale to fit the biggest. Not so easy
since in Fortran we have
     
      CHARACTER*3 SUFINP
      PARAMETER( SUFINP = 'nec' ) 

while C needs

    char SUFINP[4];
    SUFINP="nec";

Could do this with a built source but that would really mess up
portability to none Unix/POSIX systems. Could define all suffices in C
and pass to Fortran.  In some ways this solves many problems but adds
to the already massive argument lists. Could call a F77 routine from
`cnec2d` which returns the biggest length of the suffices in
`nec2d.inc`. Maybe this is the best way to do it! A second call could
then get the actual name of the input file suffix.

      INTEGER max_suf_len;
      char *buf;
      getlen_( &max_suf_len );
      buf = (char*)malloc( max_suf_len + 1 );
      getinp_( buf, max_suf_len );

### Is there a better way of doing the dynamic memory allocation?

Could use compiler specific Fortran extensions or see if possible to
call C routines from within F77 to do the allocations (this would
realistically need to work without make changes to the way the arrays
are used).

### Alter cnec2d interface to allow memory to be left allocated after call 

Use static global pointers and reused or, possible with a remalloc
in a subsequent call. Also want a call to just free existing
memory. Need a integer control param `CNEC2D_REALLOC_FREE`,
`CNED2D_REALLOC_HOLD`, `CNEC2D_ALLOC_FREE`, `CNEC2D_ALLOC_HOLD`,
`CNEC2D_FREE`.

### Using the portable cfortran interface would be good except:

1. License is restrictive. The next version is meant to be GPL'd.

2. I need about > 90 dummy args which is going to blow up a lot of brain-dead
   C preprocessors. Don't know if it is worth even beginning to consider
   working around. 

Next GPL'd version may be feasible if number of args up to 100 is possible.
 
## Allow files with ^M^J (i.e. windows format) to work on unix

And visa-versa. Not sure how easy this is in Fortran.

### Add code to check the condition number of the impedance matrix

Warn if there is likely to be a problem. See the file
`condition_number.txt` in maintainer/doc.

### Port completely to C99

Some versions already out there have done this. 

### Port to F90/95 to get dynamic arrays directly in fortran
 

## Done

 * Release objectives for 0.8.0:

   1. Initial port to cmake build system.

   2. Basic working build, install and run on Linux.

   3. Updated documentation and copyright.

   4. Tested on Debian GNU/Linux.

   5. Initial commit to mercurial repository.

   6. Review license.

 * How to force specific blas libraries

   Maybe use options:

   FIND_BLAS  = ON/OFF
   FIND_LAPACK = ON/OFF

   Default ON, as now. Off stops find_package and relies
   on user speficying via LDFLAGS, CFLAGS, FFLAGS.

   LDFLAGS=" -lblas " cmake -DSOLVER=BLAS -DFIND_BLAS=NO ../cnec2/ ; make

 * Complete overhaul of ReadMe.txt.

 * Test on Linux Intel/Athlon with and without BLAS.

 * GPL v2 or v3 or dual license.

   lapack source files.

 * Binary package:

   mlapack used: GPL3 + LAPACK 
   
   Just GPL3 if not used.

 * Update all licence headers to GPL3!

 * Review support in acinclude.m4 for auto-detecting BLAS and LAPACK
   options. Need to consider:

   + Debian sarge atlas: -llapack_atlas
   + New intel mkl: -lmkl_lapack -lmkl_ia32 -lguide -lpthread

 * Update copyright notices.

 * Add CALL FLUSH(CHIMPD) to impedance file printing code. 
   Added calls but commented out at moment.

   Need to check for FLUSH(INT) intrinsic in autoconf 
   and make netwk.f into a preprocessed file. 
 
 * Update automake/autoconf support:

   - automake 1.9.5

   - autoconf 2.59

   - solve all warnings from automake/autoconf.

     AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])], [ACTION-IF-TRUE], [ACTION-IF-FALSE])

 * Review whether I should be using @FLIBS@ from F77_LIBRARY_LDFLAGS
   when linking or is this only for C++.

 * See if the new AC_F77_MAIN and AC_F77_DUMMY_MAIN macros are
   useful and can replace my NOMAIN macro. Need autoconf 2.52.
   
   Possibly useful but don't remove need for main supression options
   with the Intel Fortran 8.1 compiler at least.

 * Devise another test or two for the built in somnec support.
   Maybe the 900 MHZ SAR validation model S11.

 * Somnec is busted on Debian sarge with gcc 3.3.5
   (same source OK on woody!)

   - Problem with bessel/hankel functions routines? Could replace them
     with routines from slatec but need quite a lot of them. May be best to 
     trace the problem. The routines shipped with NEC have some issues which
     I thought I had fixed once, but obviously not.

   - Check which variables are saved 
     in intrp in nec2++. Done. appears to fix it.

   - Also look at rom2, hankel/bessel rountines. 

   - Also do same in nec2c re incorporation of
     somnec into nec and any fixes. OK. A couple of extra SAVE's added 
     commented out for now.

   - Check NEC list archive for somnec fixes.

 * Add support for intel icc/ifort and mkl:

   - problem with missing g2c.h in c files.
     misrecognition of icc as gcc (maybe fixed
     with new autoconf). 
    
     icc must declare __USING_GCC__

   - fortran main suppression: -nofor_main

     In acinclude:

     if test "$idf_f77_nomain" = "no"; then
        strng=`${F77='f77'} -version 2>&1| grep "Intel.*Fortran"`
   
        if test -n "${strng}" ; then
           idf_f77_nomain="yes"
           FORTRAN_NOMAIN="-nofor_main"
        fi
     fi

   - with MKL Lapack

     CC=icc ./configure --with-lapack="-L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread"

   - notes on using icc/ifort.

 * Fixing missing code paths for impedance file 
   printing.  Look for write in netwk.f.

 * Update copyrights to 1998-2002.

 * Consider name change to cnec2. 

 * Upddate copyrights to 1999-2001.

 * Do I really need AM_CONFIG_HEADER? YES!

 * Update to autoconf 2.52. Run autoupdate and autoscan. See if
   my configure script is up to date. Does automake 1.4-p5 support 
   configure.ac? 

* Make blas an ADDITIONAL solver?

  cnec2 --solver=SOLVER | -r=SOLVER   SOLVER = nec2lu | blaslu
 
  SOLVER = 0    nec2lu
         = 1    blaslu

  solvers[] =                
  {
     "nec2lu", 0
     "blaslu", 1
     ""
  };

   Use global variable solver in C and pass INTEGER SOLVER to
   F77. Default is blaslu if available. Probably don't want to make code
   dependent on BLAS though. Would need to make fail gracefully if
   BLAS solver not built-in. NOT DONE.

 * Check for cygwin and mingw using AC_CANONICAL_HOST and

   case $host_os in
     *mingw32* ) MINGW32=yes;;
             * ) MINGW32=no;;
   esac
   case $host_os in
     *cygwin* ) CYGWIN=yes;;
            * ) CYGWIN=no;;
   esac

   Actually do I really need to do this. Maybe bset to remove all of
   this.

 * Don't install somnec manpage if binary not built.

 * Added checks when using built in somnec to see if ground
   arrays need to be recalculated. Update if SIG, EPSR or FMHZ
   has changed since last call.

 * Fully integrate somnec as subroutine called by NEC when necessary. This
   would be very convenient but coould be a significant undertaking. Someone
   on the NEC mailing list claims to have done this. UPDATE: Actually looks
   fairly straightforward to do this is a transparent and general way. 
   Initially lets build a libsomnec.a library in the somnec directory and have
   a configure option --enable-builtin-somnec which causes the somnec code
   to be linked into nec. Otherwise build it externally as usual. Need to 
   reinstate the ggrid.inc in the somnec.f file! The files timer.F and
   test.f are defined by both somnec and nec but they are the same. How do
   I stop somnec executable been built when built-in somnec is enabled?
   See automake manual pp.16-17 and GNU cpio example.

   So if external somnec
       build library containing all sources except test.f and timer.F
       create somnec using C/Fortran interface, test.f timer.F and library
       don't link cnec2 against libsomnec
       install somnec man page
   else
       build library containing all sources except test.f and timer.F
       don't create somnec binaries
       link nec against libsomnec
       don't install somnec man page

   Also need to modify testsuite so doesn't try to create somnec file
   when it is builtin.

   Initially default is to build external somnec until I'm happy it
   works.

 * Maybe move gnu to lib and testsuite to tests to be more standards
   compliant.

 * Put []'s into help strings for --with options when I've worked out
   what is going wrong. UPDATE: Known issue with the autoconf
   AC_HELP_STRING code. changequote(<<,>>)...changequote([,]) doesn't
   work either. Wait for upstream fix.

 * Consider integrating with my prenec program as a driver for doing
   pre-processing, running nec and post processing. In view of copyrights 
   (prenec is GPL, cnec2 is unknown) probably best to do this via a
   fork/wait interface. NOT DONE.

 * Improve accuracy and/or speed (compile time options) of the
   integration routines - may be worthwhile now other precision problems
   are fixed. Stuart says Gaussian quadrature with a low number of points
   is adequate in most cases. There a specialised integrators for
   oscillatory functions in the slatec library

   int a-->b of f(x) cos(w x)
   int a-->b of f(x) sin(w x)

   The integrations can be put into this form though there are
   complications as it involves a many to one mapping of the z plane onto
   the r plane. Also the integrand for both real and imaginary parts
   becomes singular for self-interaction terms. Actually if model is
   properly segmented you will not see the oscillations in the integrand
   over the segment length so these routines are unlikely to be of
   benefit.

   UPDATE: Analysis of INTX on a typical model with 1/10th wavelength segments
   suggests that GF is typically called 3 or 5 times per integration -
   this is the 3 and 5 point Romberg bits of INTX. GF is rarely called
   more than 5 times except for the self terms. New version of INTX using
   4 point Gaussian-Legendre Quadrature (A+S p916) was implemented. The
   improvement in speed on most of my benchmarks is only just above the
   noise floor. Suggests not worth it for well segmented models. For
   models using 1/5 wavelength segmentation or worse there may be a more
   significant difference, but the Gaussian quadrature will produce less
   accurate results with uncontrolled errors in an already borderline
   model.

   All this suggests that the number of calls to GF is already fairly
   optimal. To improve the performance of the fill without reducing
   accuracy we need to look at ways of speeding up the calculation of GF
   itself. May be possible to use splines to optimise the calculation but
   the parameter space is huge.

   There is a possible optimisation involving replacing two divides by
   a single divide, two multiplications and an assignment. Initial tests
   suggest that this makes little difference on Intel hardware
  (Improvement of about 1.7% in fill time). Divide would have to be
   really poor for this to have a significant effect. NOT DONE.

 * Investigate possibility of supplying code to put Intel Pentium FPUs
   into 53-bit IEEE double precision compliant mode. Use Bill
   Metzenthen's (http://www.linuxsupportline.com/~billm) C9X fenv
   implementation or the implementation in glib 2.1. Use autoconf to
   check for fenv.h and add appropriate code to nec.c. This may require
   linking with another library if not supported directly in the C
   library. NOT DONE.

 * Added TOTAL to column heading of total field calculations in
   nfpat. This will affect the test suite. Any tests which use NE/NH
   cards will have to have their validation files regenerated and
   CHECKED. This may not be worth the hassel. NOT DONE.

 * Update documentation for somnec.

 * Test suite script will need to be updated for new somnec.

 * Tidy up somnec - C front end, full precision constants where
   possible. Can use autoconf support in this case. Make all
   parameters specifiable on command line and allow silent operation.
   
   somnec [options] jobname[.som]

   --help           -h
   --version        -v
   --verbose        -V
   --debug=NUM      -d NUM
   --epsilon=X      -e X
   --sigma=X        -s X    X=floating point number
   --frequency=X    -f X

   Need to use strtod.

 * Bring cnec2 man page up to date. Mainly list basic usage, licence,
   and sources of further information (there are some new options I
   think).  The antennax web articles are useful.

 * Improve manual pages with links to nec manuals on the web and notes
   on the interface compatibility. Consider removing references to
   non-free software from the README and also to non-free
   documentation. What is the status of the on-line NEC manuals?

 * Add support to autoconf to test for Intel Math Kernel Library
   support. Add documentation on how to use the Intel MKL library.

   -L/opt/intel/mkl/LIB -lmkl32_lapack -lmkl32_p3 -lpthread

 * Lets not ship mkxxx scripts in nec/somnec directories.

 * Add C code to remove any existing output files before calling the
   Fortran code. Need to do this as portably as possible. POSIX. Windows? 
   Remember the NGF files - these may be read or written!

 * Improve file name suffix detection in fdrive.f to cope with
   nec.nec.nec for example.
 
 * Check blanks/zero in admittance fields. Would it be a good idea to
   explicitly zero the arrays in parsit at the start of the routine?
   Probably.

 * Audit runtests.in and acinclude.m4 for non-portable shell code. In
   particular quoting of variables, etc. Good guide in new autoconf
   manual. How portable is use of cp -f.

 * Work around need for strstr or use a replacement from GNU C
   library.

 * When checking for lapack, if we fail and BLAS_LDADD is not set
   suggest using --with-blas as well in error message.

 * Messages checking for blas/lapack at end of macros should really be 
   checking for working blas/lapack.

 * Security audit the C code for buffer overflows. Use xmalloc and
   other GNU functions where appropriate.

 * I need to actually use the function and header checks to protect my
   headers. At the moment the checks won't help much on a non-ANSI C
   system. See autoconf manual. Need to give careful consideration to
   how the gnu library is used - most of it is interlinked. See
   autoconf manual p30, 41. GNU testutils looks like a reasonable
   template to use to get the header stuff and LIBOBJS files
   correct. Sort out gnu directory, maybe also renaming lib, to make
   sure I'm providing what I need to and nothing more.

 * Check that my usage of macros like AC_CHECK_LIB is correct. It builds
   LIBS left to right - but only when no default action if true is given. 
   Also should I be using LDADD or LDFLAGS for linking in external lapack 
   libraries? Tricky one since if configure decides to use the internal
   lapack then LDADD is preferred over LDFLAGS, otherwise its the other
   way around. Could be fixed, but is it worth it?

 * When I understand how it is now supposed to work use the CYGWIN and
   MINGW32 variables in my IDF_MINGW macro. Originally set by AC_CYGWIN, 
   then went to AC_EXEEXT but this is also now obsolete. I think they are 
   now set by AC_PROG_CC and related macros. YES this is correct, I can just
   use these.

 * Go thru acinclude and configure.in and improve layout and quoting. 
   There are some guidelines in the autoconf manual.

 * Rewrite IDF_BLAS and IDF_LAPACK macros to use AC_CHECK_LIBS to
   search for blas and lapack libraries when they are requested but not
   specified:
                   AC_F77_FUNC([ZGEMM])
                   AC_F77_FUNC([ZGETRF])

                   AC_CHECK_LIB( [myblas], $ZGEMM, BLAS_LIBS="-lmyblas", , $FLIBS )

   Use a construct like that in octave. Look for atlas first. 

 * Unfortunately the autoconf 2.49b AC_F77_WRAPPER stuff won't work
   with 104 macro arguments! Thus can't do this

   #ifdef F77_FUNC
   #define NEC2D_F77( nec2d, NEC2D )
   #endif

   void NEC2D_F77( int*, .... );

   NEC2D_F77( &i,...);

   which is a bummer. Can I modify it to at least work out what the
   current name mangling is. Probably use 

   AC_F77_FUNC( NEC2D, $NEC2D )

   and then 

   AC_DEFINE( NEC2D, $NEC2D, [Mangled F77 symbol for NEC2D subroutine]).

   but may need to be careful calling the name NEC2D - maybe NEC2D_F77!

 * Use the autoconf macros for echo to sort out runtests. Replaced all
   echo's which don't need to have \n's in with multiple lines of straight
   echo. That leaves the ones which need to miss of the newline which are
   done by:

   echo $ECHO_N "Trying to do so and so...$ECHO_C"
   ....
   code
   ....
   echo "$ECHO_T done"

 * Oh grief - runtests appear to be having echo problems yet again. 
   I think we need to test how echo works in relation escape characters.

 * Add notes to README on new procedure for lapack/blas. 

 * Need to add @FLIBS@ to _LDFLAGS of hybrid C/F77 binaries.

 * Note in README on the potential limitations of the test suite -
   needs at least IEEE compliant floating point arithmetic.

 * Don't build lapack if not required. This is messy and may be best
   to just live with it. It is guaranteed to build with standard F77.
   Also need a --with-lapack=[internal|flags] to allow automatic building
   against an ATLAS equipt LAPACK library. Do this first, so

  if --with-lapack=opt is given and opt is not no we are going to use the lapack solver.
                     
      if opt is null or "internal" 

         we build the shipped mini-lapack

         if --with-blas=bopt is given 
 
            we use bopt flags as well when linking

         else

            we assume blas is found automatically

         fi 

      else 

         we try to link using opt as linker flags

      fi

   else

      we use the built in solver

   fi 

 * Bring autoconf/automake configuration in line with new CVS versions. 
   In particular make use of new Fortran support - C prototypes and 
   F77_LIBRARY_LDFLAGS.

 * Sort out the COPYING, HISTORY and AUTHORS files. Outstanding
   license concern is inclusion of LAPACK code in otherwise GPL'd
   code. Improve cnec2.1 and somnec.1 manual pages.  Reference other
   sources of information, e.g. online nec manuals for input card
   format.

 * Fix up copyright info in COPYING. Code IS GPl'ed using version 2.
   Note origin of nec2, somnec, lapack and other gpl'ed code.

 * Fix awk bug in testsuite script.

 * Add test for somnec. 20.0,0.015, 1.0,0,testsomm.som precalculate.
   Then do a RP 1 calculation on it. Actually used testvdip.nec.

 * Fix problems with somnec and g77. Localise problem: use PGI compilers
   to determine if the binary somnec file is correct when g77 uses
   -finit-local-zero. Substitute routines into g77 compiled cnec2
   one by one until source of error found.

 * Bring cnec2 man page up to date. Add brief man page for somnec.

 * Check my error messages against ones in NEC manuals. Can meaning of
   any of them be improved? Put \n in long ones! Advise on corrective
   action if possible.
g
 * May be better to add xerbla.f to lapack library to solve ASCI BLAS
   link problems on some platforms. HOWEVER will still fail blas test so
   need to include code in test macro as well. Don't know if this is a
   good idea - may be best to assume the BLAS libraries are working for
   the compiler used.
 
 * Swap meaning of -v and -V options. This affects help routine, man
   page and possibly other documentation. Add to NEWS.

 * Move all NEC2/SOMNEC code over to GPL! That just leaves the lapack
   stuff to worry about.

 * Move debug over to runtime:

   cnec2 --debug=NUM | -d=NUM    NUM=0-9 0 = non debug

   Use global variable debug in C and passed INTEGER DEBUG to F77. Don't
   use in performance critical sections! Make any files that don't need
   it .f and remove all unneeded preprocessor stuff.

 * Finish man page to extent of including all command line options, a
   brief summary of what it is and pointers to other info.
 
 * Make it so that --version --verbose dumps min, max and default array
   dimension after normal version info.

 * Remove LINEAR etc polarisation from the rdp output in separate files when
   using -m and -M.

 * testngf1.nec causes problem with awk timing calculation on some
   machines - not fill/factor time so awk processes a zero line file and
   tries to divide by NR. If NR=0 want to be zero.

 * Script runtests doesn't notice if the cnec2 executable is broken!

 * Compaq Tru64/Alpha fixes: 

   a. DEC Unix doesn't like compiling an empty file. Make lapack.F into a
      library called liblapack.a in a new directory called lapack. Remake source
      in separate files from Lapack 3.0 package. Link in with cnec2_LDADD with 
      absolute path only if blas is requested. This is cleaner anyway.

   b. Need to arrange for

      cnec2_LDFLAGS = -nofor_main

      to get defined on Tru64 Unix with DEC f77/f90/f95 compiler. Could try 
      compiling a one line Fortran file with -nofor_main and if it doesn't fail 
      assume it is required. Similarly for PGI use:

      cnec2_LDFLAGS = -Mnomain

      Config variable F77_NOMAIN.

      Algorithm:

      Try one if it works use it else try next. If none work assume don't need 
      to use an opt (e.g. gcc).

 * May be better to rename lib directory to gnu.

 * Find out how to make it so that make check finds the nec executable
   if we are not building in the source directory. May have to put tests
   in nec directory but want to avoid this. Maybe need to have a
   runtests.in with builddir subst into the path?

 * Files: If file name has a preceding directory name should the
   output go to this directory or the current directory? At the moment it
   goes to the specified directory. This affects the test suite script.

 * Benchmark version using reference BLAS from LAPACK 3.0 under
   Cygwin/mingw and Linux. Compare to builtin solver.

 * Check -m and -M do the right thing. In particular -M should
   reproduce the behaviour of Stuarts NEC.

 * Recheck array limits in NEC source manual against my source.

 * Sort out how to added autoconf tests for sizeof INTEGER and REAL*8
   (DOUBLE PRECISION). May be best to make cnec2d malloc in bytes. Then
   need to find size of Fortran INTEGER and REAL*8 types - are these
   absolutely guaranteed to be 4 and 8 bytes respectively? Don't think
   so! How do we determine the size of these?

   sizeof(REAL)             = sizeof(INTEGER)
   sizeof(COMPLEX)          = 2 * sizeof(REAL)
   sizeof(DOUBLE PRECISION) = 2 * sizeof(REAL)
             
   REAL*8 and COMPLEX*16 (DOUBLE COMPLEX) not defined. So if I can find
   sizeof(INTEGER) then I can deduce that:

   sizeof(REAL)            =     sizeof(INTEGER)
   sizeof(DOUBLE PRECSION) = 2 * sizeof(REAL)   = 2 * sizeof(INTEGER)
   sizeof(COMPLEX)         = 2 * sizeof(REAL)   = 2 * sizeof(INTEGER)

   and assume that:

   sizeof(REAL*8)     = 2 * sizeof(REAL)   = 2 * sizeof(INTEGER)
   sizeof(COMPLEX*16) = 2 * sizeof(REAL*8)   = 4 * sizeof(INTEGER)

   So how do I determine sizeof(INTEGER)?

 * Test-suite will not work if building package outside of the source
   directory. Need to look for source files in $(top_srcdir) and
   generated files in $(top_builddir). At the moment cnec2 generates
   result files in same dir as source even when doing a VPATH
   build. Simple fix is to copy input file into buiddir and run test on
   it.

 * Runtests is not portable. Remove function.

 * Write program to time *, /, sin, cos and sqrt in Fortran (no
   opt). Use results to determine any possible tuning of GF routine:

      PROGRAM TIMEOP

      INTEGER BIGNUM
      DOUBLE PRECISION A,B,C

      BIGNUM=1000000

      A=DSQRT(2.0D0)
      B=DSQRT(13.0D0)

      CALL TIMER(T1)
      DO 100 I=1,BIGNUM
100   CONTINUE
      CALL TIMER(T2)

      OVERHD = T2-T1

      CALL TIMER(T1)
      DO 200 I=1,BIGNUM
         C=A*B
200   CONTINUE
      CALL TIMER(T2)

      TIM=(T2-T1-OVERHD)/DBLE(BIGNUM)

 * See if possible to hide the autofix kludge using a wrapper script
   called configure which passes all args to configure.real and then runs
   autofix. This may upset autoconf - need to check both build in src dir
   and build in another dir. May need to maintain timestamps of script
   and real configure to make it work. May be best to have a myconf
   script and change to install instructions. Keep eye out for fixed
   version of automake.

 * Include a hard core size and number of segments/patches limit in
   the C code. Limit core to 1.6Gb, LD to 10000.

 * Diff my lapack 2.0 routines against new lapack 3.0 version for any
   improvements.

 * Rename MAXSEG to LD and check for duplicated passing to
   routines. Also change MAXSEG in c code to LD. This is more
   logical. Use MAXSEG to give PATOFF parameter.

 * Remake interface docs with gxchk and update README (info on two
   multi file output options).

 * Update version to 0.5.

 * Create binaries.
 
 * Consider adding an option to have multi file output but to same
   output files. --single-multi-file, -M. MULTIF=2.

 * Make all my Fortran an C macros use

   AC_LANG_SAVE
   AC_LANG_FORTRAN77/AC_LANG_C
   ............
   AC_LANG_RESTORE

 * Make multi file output behave like Stuarts version of NEC. This is
   short term fix. Long term move to different files for each call to
   nfpat, rdpat.

 * Make sure all files properly opened and closed.
 
 * Add documentation on compiling with cygwin/mingw for Win32 using
   reference BLAS libraries.

 * Add documentation on how to use this version compared with standard
   NEC2. Command line interface, how core size works.

 * Improve documentation in README for people not familiar with
   compiling Unix software.

 * Check for output to CHSTOT or * from nec2d.F. Remove and use IFAIL
   if necessary.

 * There is an array overflow caused by MOVE which does not check a
   priori if there is enough space in the arrays before generating new
   structures.  This is in original version but you seem to be able to
   get away with using the Fortran driver. Bombs out in C. This needs
   fixing.

 * Improve error reporting. All IFAIL=1 in nec2d.F should be replace with
   IFAIL=1,2,3,4,5 where number corresponds to error. Use an enum these
   integer to CNEC2D_ERR_XXX symbols. Write a function 

   enum {                                   IFAIL
      CNEC2D_ERR_OK,                          0
      .............                           1
      CNEC2D_ERR_RANGE,                       2
      MAXERR                                  -

   }

   static char err_msg[MAXERR][] = {
          "successful completion",
          "connection error",
   };
      
   const char *cnec2d_error( const int err_code )
   {
      if( err_code < MAXERR ) {
         return err_msg[err_code];
      }
      else {
         fprintf( stderr, "internal error: please report" );
         exit(1);
      }

   }

 * Fix bug in test suite runtests which generates ugly output if more
   than one frequency point is used. Want to average over multiple FILL
   and FACTOR values.

 * Make all STOP's a call to exit? Need to replace with a Fortran exit
   routines which can just do a stop or call exit. OR unwind all errors
   using IFAIL flag. What about performance critical routines?

 * Added INTEGER MULTIF to cnec2 args and pass to any routines with
   #ifdef MULTI_FILE. Check MULTIF not used in cnec2. Replace:
   
     #ifdef MULTI_FILE                   IF( MULTIF.EQ.1 ) THEN
     ............               with
     #endif                              ENDIF

 * C driver and front end in ./nec/

   cdriver.c
   cdriver.h

   Implement a C interface to cnec2 which pass the values of MAXSEG,
   IRESRV, NSMAX, LOADMX, NETMX, JMAX, MXCOUP, NSMAXX, NPMAX and FILNAM
   and MULTIF. Sanity check:     

   sizeof( INTEGER )    = 4
   sizeof( REAL_8 )     = 8
   sizeof( COMPLEX_16 ) = 16 

   cnec.c
   cnec.h

   Front end which uses popt to parse command line:

   Numerical Electromagnetic Code Version 2 (NEC2). Call syntax:

   cnec2 [options] jobname1[.nec] [ jobname2[.nec] ....... ]

   where jobnameN.nec are one or more nec input files to be processed. 
   The ".nec" suffix may be omitted but must be present in the actual
   file name. The valid options are:

   -h, --help            Give usage information and exit
   -v, --version         Give version information and exit
   -s, --segments=N      Maximum number of segments is N (default 300)           N
   -p, --patches=M       Maximum number of patches is M (default 0)              M
   -c, --core-size=N     Limit core size to N bytes. The suffices b, k
                         and m may be used for bytes, kilobytes and
                         megabytes respectively (default to exactly fit job)     B
   -e, --sources=N       Maximum number of sources is N (default 30)             NSMAX
   -l, --loads=N         Maximum number of loads is N ( default 30)              LOADMX
   -n, --networks=N      Maximum number of networks is N (default 10)            NETMX
   -x, --connect=N       Maximum number of segments meeting at a single  
                         wire junction is N (default 30)                         JMAX
   -u, --coupling=N      Maximum number of points in coupling
                         calculation is N (default 5)                            MXCOUP
   -S, --ngf-segments=N  Maximum number of NGF segments to which new
                         segments or patches connect is N (default 50)           NSMAXX
   -P, --ngf-segments=N  Maximum number of NGF patches to which new
                         segments or patches connect is N (default 50)           NPMAX
   -f, --frequencies=N   Maximum number of frequencies in normalised
                         impedance calculation is N (default 200)                F
   -a, --angles=N        Maximum number angles in received signal
                         strength is N (default 200))                            A
   -m, --multi-file      Enable multiple file output

   NORMF = max( F, A )
   MAXSEG = N + M

   IRESRV = ( N + 2 * M )**2    or   IRESRV = CEIL( B/16 )
                                       IRESRV >= 2 * ( N + 2 * M )
 * Redo routine interface docs with spag.

 * Remove unused variables from all routines. Need to recurse thru
   spag a few times.

 * Make main a subroutine SUBROUTINE NEC2D(.....) and create C front
   end using cfortran portable interface.

 * Dynamic memory allocation:

   a. Find out how many arguments subroutines can take in F77 and C89.

   b. Resolve all clashes between common and local variable names. Use
      PLUSFORT gxchk.out analysis.

      Go through routine by routine.

      For each variable clash

         1. Decide new name.
         2. Check for use of new name as a local in that routine or
            as a common variable.
         3. Query replace - carefully.
         4. Compile -Wall 
         5. Add ChangeLog entry. 

   c. Identify all arrays which depend on MAXSEG and IRESRV.

   d. Remove all arrays in 2 from common blocks and pass around through
      subroutine arg lists. Remember ones not in common!

 * Add --with-multifile-output to enable rdpat, near_e, near_h and
   impedance.dat files. Better names for these files
          
   impedance.dat
   near_e.dat
   near_h.dat
   rdpat.dat 

 * Add more stringent tests to check validity of lapack code with:

   a. Networks in the NEC file.
   b. A symmetrical structure.
   c. NGF.

 * Determine which LAPACK routines are actually needed and create a
   mini lapack library locally in lib dir. Make a configure option to use
   this if requested or if use-lapack given and can't find a system
   version.  BLAS remains external.

   --without-lapack

       default
       use existing code
   --with-lapack[=ldflags]
       
       use lapack routines
       define USE_LAPACK

       if flags not given 
           if -llapack -lblas works
               use this -llapack -lblas
           else 
               use builtin -L../lib -lmlapack -lblas
       else if flags=builtin 
           use builtin -L../lib -lmlapack -lblas
       else 
           use those supplied
           check works with simple program else fail.
 
   How do I prevent building local libmlapack.a in lib if not required
   using autmake.

 * Add support in nec.F to ignore blank lines and lines whose first
   character is #, ; of % in the input file.

 * Added copyright info to all my files. Not sure what to use.
   Copyright 1999, I D Flintoft <idf1@ohm.york.ac.uk>

 * Check COM save.inc - NEC3 uses CHAR array! Yes should be a
   character array (2D).

 * Trace use of files - which are opened and closed and where.  Make
   all READ/WRITE/PRINT statements used parameterised integers for the
   unit specifiers CHOT, CHIN , CHER , CHNGF ,.... which are set in
   nec.inc.
     
   CHSTIN    stdin                         5     -
   CHSTOT    stdout                        6     -
   CHSTER    stderr                        7     -
   CHINPT    input nec file                2     INFILE
   CHRSLT    output results file           3     OUTFIL
   CHSOMN    somnec file SOM2D.NEC        21     SOM2D.NEC
   CHNGFL    ngf file, IGFL               20     NGF2D.NEC
   CHTMP1    Temp file                    11     NEC11.TMP
   CHTMP2    Temp file                    12     NEC12.TMP
   CHTMP3    Temp file                    13     NEC13.TMP
   CHTMP4    Temp file                    14     NEC14.TMP
   CHTMP5    Temp file                    15     NEC15.TMP
   CHTMP6    Temp file                    16     NEC16.TMP

 * Modify Intel URL in report to http://developer.intel.com/vtune/perflibst

 * Make a simple benchmark/test suite directory. Script to run the
   generated NEC file and collect timing info and test against a known
   results file.

 * Rename SECOND to RTIME or something. There is a second intrinsic
   which we may want to use for timing if etime is not
   available. Autoconf check for etime, then second then cpu_time. Prefer
   in that order. Etime is best as it gives user time not elapse time!

   if etime
      define HAVE_ETIME
   else if second
      define HAVE_SECOND
   else if cpu_time
      define HAVE_CPU_TIME
   else
      nop

   Need to do this for nec and somnec.

 * Incorporate somnec into distribution. Need to clean up with
   pfort. Separate directory. somnec directory.

 * Autoconf checks for f77 extensions:

   + REAL*8 and COMPLEX*16 data types
   + Intrinsics: DREAL DIMAG DCMPLX DCONJG CDABS CDSQRT CDEXP DFLOAT
   + Preprocessing of .F files. Talk of adding this support to autoconf!
   + Include statement
   + List directed IO in internal files
   + Equivalencing of INTEGER and REAL*8 ok:

     INTEGER I1 , I2
     REAL*8  R1
     COMMON /FRED/ I1 ,  I2
     EQUIVALENCE (I1,R1)            

 * Put my macros into a separate m4 file. Do away with m4 directory.
     
 * Autoconf time functions. Check for working ETIME and set HAVE_ETIME
otherwise either dump time and set to zero.

 * Need to make lapack stuff autoconfable:

   --without-lapack
        default
        use existing code
   --with-lapack[=ldflags]
        use lapack routines
        define USE_LAPACK
        if flags not given assume link with -llapack -lblas
        else use those supplied
        check works with simple program else fail.

 * Linux with egcs 1.1.1 and later and maybe g77 0.5.24 support the
   ETIME intrinsic as a function or subroutine. See the LAPACK timing
   functions.

      REAL T1
      REAL TARRAY(2)
      REAL ETIME
      EXTERNAL ETIME
      T1 = ETIME(TARRAY)
      CPUSECS = TARRAY(1)

   This works on the suns too (Solaris 2.7). Make appropriate
   modifications to second.f.

 * Modify to use autoconf 2.13 and automake 1.4 with built in Fortran
   support. Ditch f2c support.

 * Look into replacing LU decomposition and back substitution with
   LAPACK/BLAS library routines - this allows PPro/PII and SMP optimised
   libraries to be used on Linux/x86 and other platforms. There is some
   info on the NEC mailing list archive on this.

	From - Mon Mar 29 08:20:35 1999
	Return-path: <davem@fs1.ece.ubc.ca>
	Envelope-to: idf1@ohm.york.ac.uk
	Delivery-date: Thu, 25 Mar 1999 22:28:30 +0000
	Received: from fs3.ece.ubc.ca ([137.82.52.241])
		by glenlivet.ohm.york.ac.uk with esmtp (Exim 2.10 #1)
		id 10QIc1-0000O8-00
		for idf1@ohm.york.ac.uk; Thu, 25 Mar 1999 22:28:29 +0000
	Received: from fs1.ece.ubc.ca (fs1.ece.ubc.ca [137.82.52.2])
		by fs3.ece.ubc.ca (8.9.0/8.9.0) with SMTP id MAA14702
		for <nec@mail.ee.ubc.ca>; Thu, 25 Mar 1999 12:29:46 -0800 (PST)
	Received: by fs1.ece.ubc.ca (4.1/SMI-4.0)
		id AA04659; Thu, 25 Mar 99 12:30:36 PST
	Message-Id: <9903252030.AA04659@fs1.ece.ubc.ca>
	Date: Thu, 25 Mar 1999 16:51:52 +0100
	From: "Juergen v.Hagen" <vonhagen@ihefiji.etec.uni-karlsruhe.de>
	To: nec-list@ee.ubc.ca
	Subject: NEC-LIST: LAPACK LU was Re: NEC-LIST: Xeon vs Pentium II - 450 
	 comparison???
	Sender: davem@ece.ubc.ca
	Content-Length: 2103

	> If performance of NEC is what you want, replace piece of Fortran code
	> which does the matrix factorization with calls to an optimized math-
	> library (e.g. Lapack). With the (free download) Intel version you get
	> a performance boost of about 3.
	> 
	> The difference between Linux and MsW is much smaller, although like
	> Dave, I would opt for Linux (and then use Greg Henry's math library.)
	> 
	> Actually, maybe someone already has rewritten the (small) piece of
	> factorization code for the speed improvement. Please let us know!

	jep. Basically you need to replace one call of FACTRS with ZGETRF, and
	one call of SOLVE with ZGETRS. Should be pretty forward. I observed
	(cf. my paper at ACES: for 1944 unknowns 378s original, 72.76 Lapack's
	LU) a speedup of more than 5 on AIX/IBM SP with the ESSL library. On a
	DEC Alpha with optimized ZGEMM speed up was at least 3. I don't use
	INTELs so I can't tell what happens there. I would be very interested
	though to hear about it.

	Interestingly, ZGETRS is able to solve for a TRANSPOSED matrix, so you
	don't even have to transpose the matrix before factorizing it.  And as
	the factorizing is the same regardless of Z^T or Z, you save the time
	for the transpose. For more info about ZGETRF look at the Users' guide
	to LAPACK and the various Lapack Working Notices available at
	www.netlib.org

	cheers
	juergen

   Need to change call to FACTR in FACTRS to

      EXTERNAL ZGETRF
      INTEGER INFO
      CALL ZGETRF(NP,NP,A(KA),NROW,IP(KA),INFO)

   and SOLVE in SOLVES to

      EXTERNAL ZGETRS
      INTEGER INFO
      CALL ZGETRS('T',NPEQ,1,A(IB),NROW,IP(IA),B(IA,IC),NROW,INFO)
	
   where INFO is a local INTEGER type. Remember external statements.

   YES - tested and it works - but need to see if worth it.

 * Reorganise - put nec in nec directory.


NEC2 Tidy And Fix Table
=======================


--------------------------------------------------------------------
              Cleaned   XREF NEC81  FTNCHEK  COMPILE  XREF NEC81-IDF
--------------------------------------------------------------------
arc.f            X         X           X       X         X
atgn2.f          X         X           X       X         X
blckin.f         X         X           X       X         X
blckot.f         X         X           X       X         X
cabc.f           X         X           X       X         X
cang.f           X         X           X       X         X
cmngf.f          X         X           X       X         X
cmset.f          X         X           X       X         X
cmss.f           X         X           X       X         X
cmsw.f           X         X           X       X         X
cmws.f           X         X           X       X         X
cmww.f           X         X           X       X         X
conect.f         X         X           X       X         X
couple.f         X         X           X       X         X
datagn.f         X         X           X       X         X
db10.f           X         X           X       X         X
db20.f           X         X           X       X         X
efld.f           X         X           X       X         X
eksc.f           X         X           X       X         X
ekscx.f          X         X           X       X         X
error.f          X         X           X       X         X
etmns.f          X         X           X       X         X
facgf.f          X         X           X       X         X
facio.f          X         X           X       X         X
factr.f          X         X           X       X         X
factrs.f         X         X           X       X         X
fbar.f           X         X           X       X         X
fblock.f         X         X           X       X         X
fbngf.f          X         X           X       X         X
ffld.f           X         X           X       X         X
fflds.f          X         X           X       X         X
gf.f             X         X           X       X         X
gfil.f           X         X           X       X         X
gfld.f           X         X           X       X         X
gfout.f          X         X           X       X         X
gh.f             X         X           X       X         X
gwave.f          X         X           X       X         X
gx.f             X         X           X       X         X
move.f           X         X           X       X         X
nefld.f          X         X           X       X         X
netwk.f          X         X           X       X         X
nfpat.f          X         X           X       X         X
nhfld.f          X         X           X       X         X
parsit.f         X         X           X       X         X
patch.f          X         X           X       X         X
pcint.f          X         X           X       X         X
prnt.f           X         X           X       X         X
qdsrc.f          X         X           X       X         X
rdpat.f          X         X           X       X         X
readgm.f         X         X           X       X         X
readmn.f         X         X           X       X         X
reblk.f          X         X           X       X         X
reflc.f          X         X           X       X         X
rom2.f           X         X           X       X         X
sbf.f            X         X           X       X         X
second.f         X         X           X       X         X
sflds.f          X         X           X       X         X
solgf.f          X         X           X       X         X
solve.f          X         X           X       X         X
solves.f         X         X           X       X         X
subph.f          X         X           X       X         X 
tbf.f            X         X           X       X         X
test.f           X         X           X       X         X
trio.f           X         X           X       X         X
unere.f          X         X           X       X         X
upcase.f         X         X           X       X         X
wire.f           X         X           X       X         X
zint.f           X         X           X       X         X
----------------------------------------------------------------------------




