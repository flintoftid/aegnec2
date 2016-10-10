C* 
C* aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
C* Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>.
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
 
      SUBROUTINE GX( ZZ , RH , XK , GZ , GZP )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Segment end contributions for thin wire approximation kernel.      **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ZZ                                                        **
C* INPUT  - RH                                                        **
C* INPUT  - XK                                                        **
C* OUTPUT - GZ                                                        **
C* OUTPUT - GZP                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   EKSC    EKSCX                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 ZZ , RH , XK
      COMPLEX*16 GZ , GZP

C     Local variables.
      REAL*8 R , RR2 , RK
 

      RR2 = ZZ*ZZ + RH*RH
      R = DSQRT(RR2)
      RK = XK*R
      GZ = DCMPLX(DCOS(RK),-DSIN(RK))/R
      GZP = -DCMPLX(1.0D0,RK)*GZ/RR2
 
      RETURN
 
      END
