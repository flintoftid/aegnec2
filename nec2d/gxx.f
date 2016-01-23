C* 
C* aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
C* Copyright (C) 1998-2016 Ian David Flintoft <ian.flintoft@york.ac.uk>
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

      SUBROUTINE GXX( ZZ , RH , A , A2 , XK , IRA , G1 , G1P , G2 , 
     &                G2P , G3 ,GZP )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Segment end contributions for extended thin wire approximation     **
C* kernel.                                                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ZZ                                                        **
C* INPUT  - RH                                                        **
C* INPUT  - A                                                         **
C* INPUT  - A2                                                        **
C* INPUT  - XK                                                        **
C* INPUT  - IRA                                                       **
C* OUTPUT - G1                                                        **
C* OUTPUT - G1P                                                       **
C* OUTPUT - G2                                                        **
C* OUTPUT - G2P                                                       **
C* OUTPUT - G3                                                        **
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
C* called by   EKSCX                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER IRA
      REAL*8 ZZ , RH , A , A2 , XK
      COMPLEX*16 G1 , G1P , G2 , G2P , G3 , GZP

C     Local variables.
      REAL*8 R , RR2 , R4 , RH2 , RK , RK2 , TT1 , TT2
      COMPLEX*16 GZ , C1 , C2 , C3


      RR2 = ZZ*ZZ + RH*RH
      R = DSQRT(RR2)
      R4 = RR2*RR2
      RK = XK*R
      RK2 = RK*RK
      RH2 = RH*RH
      TT1 = 0.25D0*A2*RH2/R4
      TT2 = 0.5D0*A2/RR2
      C1 = DCMPLX(1.0D0,RK)
      C2 = 3.0D0*C1 - RK2
      C3 = DCMPLX(6.0D0,RK)*RK2 - 15.0D0*C1
      GZ = DCMPLX(DCOS(RK),-DSIN(RK))/R
      G2 = GZ*(1.0D0+TT1*C2)
      G1 = G2 - TT2*C1*GZ
      GZ = GZ/RR2
      G2P = GZ*(TT1*C3-C1)
      GZP = TT2*C2*GZ
      G3 = G2P + GZP
      G1P = G3*ZZ
      IF ( IRA.EQ.1 ) GOTO 2
      G3 = (G3+GZP)*RH
      GZP = -ZZ*C1*GZ
      IF ( RH.GT.1.D-10 ) GOTO 1
      G2 = (0.0D0,0.0D0)
      G2P = (0.0D0,0.0D0)
      RETURN
    1 G2 = G2/RH
      G2P = G2P*ZZ/RH
      RETURN
    2 TT2 = 0.5D0*A
      G2 = -TT2*C1*GZ
      G2P = TT2*GZ*C2/RR2
      G3 = RH2*G2P - A*GZ*C1
      G2P = G2P*ZZ
      GZP = -ZZ*C1*GZ
 
      RETURN
 
      END
