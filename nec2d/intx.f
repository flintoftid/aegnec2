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

      SUBROUTINE INTX( EL1 , EL2 , BB , IJ , SGR , SGI )

C*--------------------------------------------------------------------**
C*                                                                    **
C* INTX performs numerical integration of exp(jkr)/r by the method of **
C* variable interval width romberg integration. The integrand value   **
C* is supplied by subroutine GF.                                      **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - EL1                                                       **
C* INPUT  - EL2                                                       **
C* INPUT  - BB                                                        **
C* INPUT  - IJ                                                        **
C* OUTPUT - SGR                                                       **
C* OUTPUT - SGI                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       GF      TEST                                           **
C* called by   EKSC    EKSCX                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL GF , TEST

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER IJ
      REAL*8 EL1 , EL2 , BB , SGR , SGI

C     Local variables.
      INTEGER NM , NS , NT , NTS , NX
      REAL*8 DZ , DZOT , EP , FNM , FNS , G1I , G1R , G2I , G2R , 
     &       G3I , G3R , G4I , G4R , G5I , G5R , RX , SSS , T00I , 
     &       T00R , T01I , T01R , T02I , T02R , T10I , T10R , T11I , 
     &       T11R , T20I , T20R , TE1I , TE1R , TE2I , TE2R , ZZ , 
     &       ZE , ZEND , ZP

C     Data initialisation.
      DATA NX /1/
      DATA NM /65536/
      DATA NTS /4/
      DATA RX /1.0D-4/


      ZZ = EL1
      ZE = EL2
      IF ( IJ.EQ.0 ) ZE = 0.0D0
      SSS = ZE - ZZ
      FNM = NM
      EP = SSS/(10.0D0*FNM)
      ZEND = ZE - EP
      SGR = 0.0D0
      SGI = 0.0D0
      NS = NX
      NT = 0
      CALL GF(ZZ,G1R,G1I)
    1 FNS = NS
      DZ = SSS/FNS
      ZP = ZZ + DZ
      IF ( ZP-ZE ) 3 , 3 , 2
    2 DZ = ZE - ZZ
      IF ( DABS(DZ)-EP ) 17 , 17 , 3
    3 DZOT = DZ*0.5D0
      ZP = ZZ + DZOT
      CALL GF(ZP,G3R,G3I)
      ZP = ZZ + DZ
      CALL GF(ZP,G5R,G5I)
    4 T00R = (G1R+G5R)*DZOT
      T00I = (G1I+G5I)*DZOT
      T01R = (T00R+DZ*G3R)*0.5D0
      T01I = (T00I+DZ*G3I)*0.5D0
      T10R = (4.0D0*T01R-T00R)/3.0D0
      T10I = (4.0D0*T01I-T00I)/3.0D0

C     Test convergence of 3 point romberg result.

      CALL TEST(T01R,T10R,TE1R,T01I,T10I,TE1I,0.0D0)
      IF ( TE1I-RX ) 5 , 5 , 6
    5 IF ( TE1R-RX ) 8 , 8 , 6
    6 ZP = ZZ + DZ*0.25D0
      CALL GF(ZP,G2R,G2I)
      ZP = ZZ + DZ*0.75D0
      CALL GF(ZP,G4R,G4I)
      T02R = (T01R+DZOT*(G2R+G4R))*0.5D0
      T02I = (T01I+DZOT*(G2I+G4I))*0.5D0
      T11R = (4.0D0*T02R-T01R)/3.0D0
      T11I = (4.0D0*T02I-T01I)/3.0D0
      T20R = (16.0D0*T11R-T10R)/15.0D0
      T20I = (16.0D0*T11I-T10I)/15.0D0

C     Test convergence of 5 point romberg result.

      CALL TEST(T11R,T20R,TE2R,T11I,T20I,TE2I,0.0D0)
      IF ( TE2I-RX ) 7 , 7 , 14
    7 IF ( TE2R-RX ) 9 , 9 , 14
    8 SGR = SGR + T10R
      SGI = SGI + T10I
      NT = NT + 2
      GOTO 10
    9 SGR = SGR + T20R
      SGI = SGI + T20I
      NT = NT + 1
   10 ZZ = ZZ + DZ
      IF ( ZZ-ZEND ) 11 , 17 , 17
   11 G1R = G5R
      G1I = G5I
      IF ( NT-NTS ) 1 , 12 , 12
   12 IF ( NS-NX ) 1 , 1 , 13

C     Double step size.

   13 NS = NS/2
      NT = 1
      GOTO 1
   14 NT = 0
      IF ( NS-NM ) 16 , 15 , 15
   15 WRITE (CHRSLT,20) ZZ
      GOTO 9

C     Halve step size.

   16 NS = NS*2
      FNS = NS
      DZ = SSS/FNS
      DZOT = DZ*0.5D0
      G5R = G3R
      G5I = G3I
      G3R = G2R
      G3I = G2I
      GOTO 4
   17 CONTINUE
      IF ( IJ ) 19 , 18 , 19

C     Add contribution of near singularity for diagonal term.

   18 SGR = 2.0D0*(SGR+DLOG((DSQRT(BB*BB+SSS*SSS)+SSS)/BB))
      SGI = 2.0D0*SGI
   19 CONTINUE

      RETURN
 
   20 FORMAT (' STEP SIZE LIMITED AT Z=',F10.5)
 
      END
