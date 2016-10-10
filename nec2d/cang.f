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
C*------------------------- Dummy Arguments --------------------------**
 
      REAL*8 FUNCTION CANG( ZZ )

C*--------------------------------------------------------------------**
C*                                                                    **  
C* CANG returns the phase angle of a complex number in degrees.       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ZZ                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ATGN2                                                  **
C* called by   NEC2D   NFPAT   RDPAT                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines
      REAL*8 ATGN2
      EXTERNAL ATGN2

C     Dummy arguments.
      COMPLEX*16 ZZ


      CANG = ATGN2(DIMAG(ZZ),DBLE(ZZ))*57.29577951308232087679815D0
 
      RETURN
 
      END
