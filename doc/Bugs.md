
# Known Bugs In AEGNEC2

Please report bugs using the bitbucket issue tracker at
<https://bitbucket.org/uoyaeg/aegnec2/issues> or by email 
to <ian.flintoft@googlemail.com>.

## SOMNEC appears to be broken on 64-bit systems

This seems to be related to the C/Fortran interface.

## The validation of testngf1.nec and testngf2.nec is known to fail on alpha platforms 

Occurs with both Tru64 Unix and Linux. This is due to a precision problems in 
regenerating the geometry data from the binary NGF file causing tan(0,x) to be 
180 instead of 0 due to noise in small values of x. atgn2.f wraps the standard 
DATAN2(Y,X) function such that it returns 0 if X==Y==0. If Y==0 but X is a very 
small number then rounding in the NGF code can flip its sign causing ATGN2 to 
flip from 0 to 180. These values are effectively degenerate so the failure is 
not significant. It is however annoying and will be fixed eventually.

Workaround: Ignore errors in testsuite caused by this problem.
 
## VSRC array limit in NETWK

This limit was 30 and is now changed to NETMX. The dimension is not mentioned in 
manual section on array limits. From the code it needs to be as large as NTSC, 
the number of network voltage sources. I believe this means that the appropriate 
dimension is NETMX but keep this in mind if any network relate problems emerge.

## PATOFF offset in ICON arrays limits maximum number of segments

Originally hard coded as 10000 (9999, 9998 as well) but now set via a parameter. 
Are there any other builtin limits on the number of segments related to PATOFF 
that I've missed?


# Fixed Issues

## Strange filename buffer problem.

Dormant.

On some occasion the results appear to depend on the name of the input file. It 
appears to only occur with gcc on Linux systems and may be related to the 
C-Fortran interface. May be a bug in the nec code, the compiler or the C 
libraries. UPDATE: Only seems to happen with ASCI Red BLAS libraries, it is OK 
with altas and built-in solver. Is it a buffer problem in the ASCI Red 
libraries? Also doesn't happen with PGI compilers on Linux - is it an interface 
problem with the ASCI Red libraries? Almost certainly an ASCI RED problem. 

## The somnec program is broken when compiled with gcc-3.0. 

Fixed in 0.7.2.
 
The Bessel functions and Hankel functions routines are the source of the 
problem. This will either be fixed or worked around in a later version, possible 
by replacing the Bessel and Hankel functions routines which have proved 
problematic in the past.

## VPATH builds will not work SunOS

Fixed in 0.7.1.

Due to a portability problem in the new version of autoconf. The compile command 
comes out as

    cc -DHAVE_CONFIG_H -I. -I../../AEGNEC2-0.6.6a1/lib -I.. -I.. -I../../AEGNEC2-0.6.6a1/lib \
      `test -f ../../AEGNEC2-0.6.6a1/lib/error.c || echo '../../AEGNEC2-0.6.6a1/lib/'`error.c

instead of 

    cc -DHAVE_CONFIG_H -I. -I../../AEGNEC2-0.6.6a1/lib -I.. -I.. -I../../AEGNEC2-0.6.6a1/lib \
      `test -f ../../AEGNEC2-0.6.6a1/lib/error.c || echo '../../AEGNEC2-0.6.6a1/lib/'`../../AEGNEC2-0.6.6a1/lib/error.c

Workaround: Build in source directory.

## Post-processing of the results in the testsuite generates error messages for some of the NGF tests

Fixed in 0.6.4.

Specifically those tests without a FILL and FACTOR time in the results file) 
with some versions of awk fail. This does not actually cause real problems, but 
does look bad. Need to find a lower common denominator of awk functionality.

## On Solaris x86 and Linux using gcc 2.95 solutions using somnec are completely broken.

Fixed in 0.6.3.

OK on Linux with PGI. Part of this is caused by uninitialised variables in both 
somnec and nec (intrp.f). Some of these problems have been fixed. However, even 
compiling with -finit-local-zero on Linux does not fix the problem - it does 
stop the segmentation faults but the results are garbage.

UPDATE: hankel.f and bessel.f in somnec and intrp.f in nec require that some 
local variables are saved between calls to the subroutine. By default g77 (on 
all platforms) uses automatic local variables and so these routines fail. This 
is fixed as of 0.6.3 by issuing a SAVE in these routines to save all local 
variables. Work around in earlier versions by specifying the -fno-automatic 
option for g77.
