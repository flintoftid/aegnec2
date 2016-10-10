C* 
C* aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
C* Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>
C*
C* This program is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C*
C* This program is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C* GNU General Public License for more details.
C*
C* You should have received a copy of the GNU General Public License
C* along with this program.  If not, see <http://www.gnu.org/licenses/>.
C* 
C*--------------------------------------------------------------------**

      SUBROUTINE LAMBDA ( T , XLAM , DXLAM )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute integration parameter XLAM=LAMBDA from parameter T.        **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - T                                                         **
C* OUTPUT - XLAM                                                      **
C* OUTPUT - DXLAM                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  A       B                                              **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   ROM1    SAOA                                           **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 T
      COMPLEX*16 XLAM , DXLAM

C     Common storage.
      INCLUDE 'cntour.inc'


      DXLAM=B-A
      XLAM=A+DXLAM*T

      RETURN

      END
