
# AEG NEC2 Licence

aegnec2 - Applied Electromagnetics Group dynamically allocated 
Numerical Electromagnetics Code Version 2

Copyright (C) 1998-2016 Ian David Flintoft <ian.flintoft@googlemail.com>

This version of NEC2 is distributed under the terms of the GNU General
Public License Version 3. 

## aegnec2 source code

The code in the nec2d directory is derived from the code in
nec2_src.tar.Z on The Unofficial Numerical Electromagnetics Code (NEC)
Archives that used to be at http://www.qsl.net/wb6tpu/swindex.html. 
The SOMNEC code is derived from somnec.tar.Z at the same location. 
These codes are in the public domain. The copyright on the modified
and new code is held by Ian David Flintoft <ian.flintoft@googlemail.com>
and it is distributred under the terms of the GNU General Public 
License Version 3. See the file doc/gpl-3.0.txt in the souce code 
distribution.

## LAPACK source code

The mlapack directory contains a few routines from the 
[LAPACK 3.2.1 mathematics libraries](http://www.netlib.org/lapack) which are
built into a mini LAPACK library under certain build options. This is a 
convenience so users are not required to build a full version of LAPACK 
in order to use the BLAS based version of cnec2. These files are 
distributed under their own GPL compatible license - see the file
mlapack/[LAPACK_LICENSE.txt][] in the souce code distribution.
 

[LAPACK_LICENSE.txt]: https://github.com/flintoftid/aegnec2/blob/master/mlapack/LAPACK_LICENSE.txt
