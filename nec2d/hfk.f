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

      SUBROUTINE HFK( EL1 , EL2 , RHK , ZPKX , SGR , SGI )

C*--------------------------------------------------------------------**
C*                                                                    **
C* HFK computes the H field of a uniform current filament by          **
C* numerical integration.                                             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - EL1                                                       **
C* INPUT  - EL2                                                       **
C* INPUT  - RHK                                                       **
C* INPUT  - ZPKX                                                      **
C* OUTPUT - SGR                                                       **
C* OUTPUT - SGI                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    RHKS    ZPK2                                           **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       GH      TEST                                           **
C* called by   HSFLX                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL GH , TEST

C     Dummy arguments.
      REAL*8 EL1 , EL2 , RHK , ZPKX , SGR , SGI

C     Local variables.
      INTEGER NM , NS , NT , NTS , NX
      REAL*8 DZ , DZOT , EP , G1I , G1R , G2I , G2R , G3I , G3R , 
     &       G4I , G4R , G5I , G5R , RX , SSS , T00I , T00R , T01I , 
     &       T01R , T02I , T02R , T10I , T10R , T11I , T11R , T20I , 
     &       T20R , TE1I , TE1R , TE2I , TE2R , ZZ , ZE , ZEND , ZP

C     Common storage.
      INCLUDE 'tmh.inc'
 
C     Data initialisation.
      DATA NX  /1/
      DATA NM /65536/
      DATA NTS /4/
      DATA RX /1.0D-4/
 

      ZPK2 = ZPKX
      RHKS = RHK*RHK
      ZZ = EL1
      ZE = EL2
      SSS = ZE - ZZ
      EP = SSS/(10.0D0*NM)
      ZEND = ZE - EP
      SGR = 0.0D0
      SGI = 0.0D0
      NS = NX
      NT = 0
      CALL GH(ZZ,G1R,G1I)
    1 DZ = SSS/NS
      ZP = ZZ + DZ
      IF ( ZP-ZE ) 3 , 3 , 2
    2 DZ = ZE - ZZ
      IF ( DABS(DZ)-EP ) 17 , 17 , 3
    3 DZOT = DZ*0.5D0
      ZP = ZZ + DZOT
      CALL GH(ZP,G3R,G3I)
      ZP = ZZ + DZ
      CALL GH(ZP,G5R,G5I)
    4 T00R = (G1R+G5R)*DZOT
      T00I = (G1I+G5I)*DZOT
      T01R = (T00R+DZ*G3R)*0.5D0
      T01I = (T00I+DZ*G3I)*0.5D0
      T10R = (4.0D0*T01R-T00R)/3.0D0
      T10I = (4.0D0*T01I-T00I)/3.0D0
      CALL TEST(T01R,T10R,TE1R,T01I,T10I,TE1I,0.0D0)
      IF ( TE1I-RX ) 5 , 5 , 6
    5 IF ( TE1R-RX ) 8 , 8 , 6
    6 ZP = ZZ + DZ*0.25D0
      CALL GH(ZP,G2R,G2I)
      ZP = ZZ + DZ*0.75D0
      CALL GH(ZP,G4R,G4I)
      T02R = (T01R+DZOT*(G2R+G4R))*0.5D0
      T02I = (T01I+DZOT*(G2I+G4I))*0.5D0
      T11R = (4.0D0*T02R-T01R)/3.0D0
      T11I = (4.0D0*T02I-T01I)/3.0D0
      T20R = (16.0D0*T11R-T10R)/15.0D0
      T20I = (16.0D0*T11I-T10I)/15.0D0
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
   13 NS = NS/2
      NT = 1
      GOTO 1
   14 NT = 0
      IF ( NS-NM ) 16 , 15 , 15
   15 WRITE (3,18) ZZ
      GOTO 9
   16 NS = NS*2
      DZ = SSS/NS
      DZOT = DZ*0.5D0
      G5R = G3R
      G5I = G3I
      G3R = G2R
      G3I = G2I
      GOTO 4
   17 CONTINUE
      SGR = SGR*RHK*0.5D0
      SGI = SGI*RHK*0.5D0
 
      RETURN
 
   18 FORMAT (' STEP SIZE LIMITED AT Z=',F10.5)
 
      END
