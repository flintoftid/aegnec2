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

      SUBROUTINE EKSCX( BBX , SSS , ZZ , RHX , XK , IJ , INX1 , INX2 ,
     &                  EEZS , ERS , EEZC , ERC , EEZK , ERK )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute E field of sine, cosine, and constant current filaments by **
C* extended thin wire approximation.                                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - BBX                                                       **
C* INPUT  - SSS                                                       **
C* INPUT  - ZZ                                                        **
C* INPUT  - RHX                                                       **
C* INPUT  - XK                                                        **
C* INPUT  - IJ                                                        **
C* INPUT  - INX1                                                      **
C* INPUT  - INX2                                                      **
C* OUTPUT - EEZS                                                      **
C* OUTPUT - ERS                                                       **
C* OUTPUT - EEZC                                                      **
C* OUTPUT - ERC                                                       **
C* OUTPUT - EEZK                                                      **
C* OUTPUT - ERK                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    IJX     RKB2    ZPK                                    **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       GX      GXX     INTX                                   **
C* called by   EFLD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL GX , GXX , INTX

C     Dummy arguments.
      INTEGER IJ , INX1 , INX2
      REAL*8 BBX , SSS , ZZ , RHX , XK
      COMPLEX*16 EEZS , ERS , EEZC , ERC , EEZK , ERK

C     Local variables.
      INTEGER IRA
      REAL*8 A2 , BB , BK , BK2 , CINT , CONX(2) , CS , RH , RHK , 
     &       SH , SHK , SINT , SS , Z1 , ZZ2
      COMPLEX*16 CON , GZ1 , GZ2 , GZP1 , GZP2 , GR1 , GR2 , GRP1 , 
     &           GRP2 , GRK1 , GRK2 , GZZ1 , GZZ2

C     Common storage.
      INCLUDE 'tmi.inc'

C     Equivalent storage.
      EQUIVALENCE (CONX,CON)
 
C     Data initialisation
      DATA CONX/0.0D0,4.771345159236942258888898D0/
 

      IF ( RHX.LT.BBX ) GOTO 1
      RH = RHX
      BB = BBX
      IRA = 0
      GOTO 2
    1 RH = BBX
      BB = RHX
      IRA = 1
    2 SH = 0.5D0*SSS
      IJX = IJ
      ZPK = XK*ZZ
      RHK = XK*RH
      RKB2 = RHK*RHK
      SHK = XK*SH
      SS = DSIN(SHK)
      CS = DCOS(SHK)
      ZZ2 = SH - ZZ
      Z1 = -(SH+ZZ)
      A2 = BB*BB
      IF ( INX1.EQ.2 ) GOTO 3
      CALL GXX(Z1,RH,BB,A2,XK,IRA,GZ1,GZP1,GR1,GRP1,GRK1,GZZ1)
      GOTO 4
    3 CALL GX(Z1,RHX,XK,GZ1,GRK1)
      GZP1 = GRK1*Z1
      GR1 = GZ1/RHX
      GRP1 = GZP1/RHX
      GRK1 = GRK1*RHX
      GZZ1 = (0.0D0,0.0D0)
    4 IF ( INX2.EQ.2 ) GOTO 5
      CALL GXX(ZZ2,RH,BB,A2,XK,IRA,GZ2,GZP2,GR2,GRP2,GRK2,GZZ2)
      GOTO 6
    5 CALL GX(ZZ2,RHX,XK,GZ2,GRK2)
      GZP2 = GRK2*ZZ2
      GR2 = GZ2/RHX
      GRP2 = GZP2/RHX
      GRK2 = GRK2*RHX
      GZZ2 = (0.0D0,0.0D0)
    6 EEZS = CON*((GZ2-GZ1)*CS*XK-(GZP2+GZP1)*SS)
      EEZC = -CON*((GZ2+GZ1)*SS*XK+(GZP2-GZP1)*CS)
      ERS = -CON*((ZZ2*GRP2+Z1*GRP1+GR2+GR1)*SS-(ZZ2*GR2-Z1*GR1)*CS*XK)
      ERC = -CON*((ZZ2*GRP2-Z1*GRP1+GR2-GR1)*CS+(ZZ2*GR2+Z1*GR1)*SS*XK)
      ERK = CON*(GRK2-GRK1)
      CALL INTX(-SHK,SHK,RHK,IJ,CINT,SINT)
      BK = BB*XK
      BK2 = BK*BK*0.25D0
      EEZK = -CON*(GZP2-GZP1+XK*XK*(1.0D0-BK2)*DCMPLX(CINT,-SINT)
     &      -BK2*(GZZ2-GZZ1))
 
      RETURN
 
      END
