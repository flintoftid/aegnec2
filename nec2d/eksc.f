C* 
C* cnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
C* Copyright (C) 1998-2009 Ian David Flintoft <idf1@ohm.york.ac.uk>
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

      SUBROUTINE EKSC( SSS , ZZ , RH , XK , IJ , EEZS , ERS , EEZC , 
     &                 ERC , EEZK , ERK )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute E field of sine, cosine, and constant current filaments by **
C* thin wire approximation.                                           **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - SSS                                                       **
C* INPUT  - ZZ                                                        **
C* INPUT  - RH                                                        **
C* INPUT  - XK                                                        **
C* INPUT  - IJ                                                        **
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
C* calls       GX      INTX                                           **
C* called by   EFLD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL GX, INTX

C     Dummy arguments.
      INTEGER IJ 
      REAL*8 SSS , ZZ , RH , XK  
      COMPLEX*16 EEZS , ERS , EEZC , ERC , EEZK , ERK

C     Local variables.
      REAL*8 CINT , CS , RHK , SH , SHK , SINT , SS , Z1 , ZZ2
      REAL*8 CONX(2)
      COMPLEX*16 CON , GZ1 , GZ2 , GP1 , GP2 , GZP1 , GZP2

C     Common storage.
      INCLUDE 'tmi.inc'

C     Equivalent storage.
      EQUIVALENCE (CONX,CON)
 
C     Data initialisation.
      DATA CONX/0.0D0,4.771345159236942258888898D0/
 
 
      IJX = IJ
      ZPK = XK*ZZ
      RHK = XK*RH
      RKB2 = RHK*RHK
      SH = 0.5D0*SSS
      SHK = XK*SH
      SS = DSIN(SHK)
      CS = DCOS(SHK)
      ZZ2 = SH - ZZ
      Z1 = -(SH+ZZ)
      CALL GX(Z1,RH,XK,GZ1,GP1)
      CALL GX(ZZ2,RH,XK,GZ2,GP2)
      GZP1 = GP1*Z1
      GZP2 = GP2*ZZ2
      EEZS = CON*((GZ2-GZ1)*CS*XK-(GZP2+GZP1)*SS)
      EEZC = -CON*((GZ2+GZ1)*SS*XK+(GZP2-GZP1)*CS)
      ERK = CON*(GP2-GP1)*RH
      CALL INTX(-SHK,SHK,RHK,IJ,CINT,SINT)
      EEZK = -CON*(GZP2-GZP1+XK*XK*DCMPLX(CINT,-SINT))
      GZP1 = GZP1*Z1
      GZP2 = GZP2*ZZ2
      IF ( RH.LT.1.D-10 ) GOTO 1
      ERS = -CON*((GZP2+GZP1+GZ2+GZ1)*SS-(ZZ2*GZ2-Z1*GZ1)*CS*XK)/RH
      ERC = -CON*((GZP2-GZP1+GZ2-GZ1)*CS+(ZZ2*GZ2+Z1*GZ1)*SS*XK)/RH
      RETURN
    1 ERS = (0.0D0,0.0D0)
      ERC = (0.0D0,0.0D0)
 
      RETURN
 
      END
