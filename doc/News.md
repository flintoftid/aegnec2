
# AEG NEC2 User Visible Changes

This file notes major user visible changes to the code. Details may be
found in the ChangeLog file.

## 15/05/2009 Version 0.8.0 Released

* Initial port to a new development system: 

  + CNEC2 is CMake now configured and built using CMake. Version 2.6.4 
    or later is required.

  + The source code is now managed using the Mercurial VCS.

* The license has been changed from version 2 to version 3 of the GNU 
  General Public License.

* The mini-LAPACK library has been update from version 2.0 to version
  3.2.1.

## 25/05/2005 Version 0.7.2 Released

* Somnec is now hopefully fully fixed with respect to saved 
  variables and initialisation.

* When multi-file output is used the impedance file is flushed
  to disk after every frequency step.

## 16/07/2002 Version 0.7.1 Released

* Renamed package from nec2d to cnec2. Also the executables change
  from nec2d to cnec2 and fnec2d becomes fnec2. This was done because
  there is already a version of NEC2 out there called nec2d!

* The gnu and testsuite directory have been renamed lib and tests
  repsectively.

* SOMNEC can now be built into the cnec2 executable using the
  --enable-builtin-somnec configure option. At a later date this will 
  become the default behaviour.

* --with-mingw configure option removed. use FFLAGS instead.

* Development now based on autoconf 2.53 and automake 1.6.2.

## 26/07/2001 Version 0.7.0 Released

* Documentation updates.

* New somnec C front end. The old interface is still available in the
  executable fsomnec.

* The 2000 unknown version of nec with Fortran front end is now
  called fnec2d (was called nec2000).

* Development now requires autoconf 2.50 and automake 2.4-p5.

## 04/05/2001 Version 0.6.8 Released

* Maintenance release.

* Fixed problems with jobname extraction with Fortran front end and
  problems with file opening on Linux systems using UNKNOWN
  status. All new files now opened using 'NEW'. This requires the old
  files are deleted by the C front end before calling the Fortran
  code.

## 12/01/2001 Version 0.6.7 Released

* Maintenance release. Updates to configuration and utility code.

* Fixed feature with filename dependent results in C - Fortran
  interface.

* Output files are now opened with status 'NEW' and the old output 
  files removed in the C front end. The Fortran front end will now 
  fail if there are existing output files.

## 13/10/2000 Version 0.6.6 Released

* Maintenance release. Updates to configuration and utility code.

## 04/10/2000 Version 0.6.5 Released

* Added support for using an external LAPACK library to the configure
  script. This is useful when using ATLAS for example.

* Now unquestionable licensed under the GPL.

* Development now requires CVS versions of automake (1.4a) and
 autoconf (2.49b).

* Bug fixes.

## 17/04/2000 Version 0.6.4 Released

* Bug fixes.

## 21/03/2000 Version 0.6.3 Released

* Testsuite case for lossy ground using somnec.

* somnec now works with g77 and other compilers which don't save
  local variables, i.e. use automatic variables by default.

* The -M and -m options now work with ground wave output from the RP
  card, i.e. RP 1 .....

## 10/10/1999 Version 0.6.2 Released
 
* Improved quality of error messages.

* Fixed truncation of output files with -M options.

* Improved man page for nec2d.

* Polarisation is not printed in rdp files to aid plotting.

## 06/09/1999 Version 0.6.1 Released

* A few bug fixes, mainly for configuration on alphas.

## 05/09/1999 Version 0.6 Released

* Now supports transparent compilation with PGI compilers (versions
  3.0 and higher) on Linux and Solaris. Older versions still require
  hand linking as described in the README file).

* Now supports Digital/Compaq compilers on Compaq Tru64 Unix.

* Now uses LAPACK 3.0 LU routines. These give a slight improvement in
  performance. The LAPACK routines are now in a separate library.

* The maintainer source code in separate files now shipped in main
  distribution.

* Autofix no longer required after configuration step.

* Supports VPATH builds.

## 20/08/1999 Version 0.5 Released

* C front end for dynamic array allocation.

* Allows embedded comments in input file.

* Runtime choice of using multiple file output using two different
 approaches.

* Improved test suite.

* Support for Cygwin/mingw Win32 compilation environment.

## 09/04/1999 Version 0.3 Released

* Working timing support using ETIME, CPU_TIME or SECOND.

* Support for using LAPACK for matrix factorisation and back
  substitution.

* Added somnec.

* Added test suite.

* Improved configuration.

## 06/12/1998 Version 0.2 Released

* Major code cleanup. 
 
* Common blocks realigned to meet standards.

* User configurable parameters isolated into nec.inc.

* Fully double precision including all constants.

## 05/10/1998 Version 0.1 Released

* Initial release based on original NEC2 source code.
