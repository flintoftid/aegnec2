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

      SUBROUTINE ROM1 ( N , SUM , NX )

C*--------------------------------------------------------------------**
C*                                                                    **
C* ROM1 integrates the 6 Sommerfeld integrals from a to b in lambda.  **
C* the method of variable interval width romberg integration is used. **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - N                                                         **
C* OUTPUT - SUM                                                       **
C* INPUT  - NX                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  A       B                                              **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       LAMBDA  SAOA    TEST                                   **
C* called by   EVLUA   GSHANK                                         **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL LAMBDA , SAOA , TEST  

C     Dummy arguments.
      INTEGER N , NX
      COMPLEX*16 SUM(6)

C     Local variables.
      INTEGER I , LSTEP , NM , NOGO , NS , NT , NTS
      REAL*8 DZ , DZOT , EP , RX , S , TI , TR , Z , ZE , ZEND
      COMPLEX*16 T00 , T02 , T11
      COMPLEX*16 G1(6) , G2(6) , G3(6) , G4(6) , G5(6) , T01(6), 
     &           T10(6), T20(6)

C     Save local variables.
      SAVE Z , ZE , S , EP , ZEND , DZ , DZOT , TR , TI 
      SAVE T00 , T11 , T02
      SAVE G1 , G2 , G3 , G4 , G5 , T01 , T10 , T20

C     Common storage.
      INCLUDE 'cntour.inc'

C     Data initialisation.
      DATA NM /131072/
      DATA NTS /4/
      DATA RX /1.0D-4/

      LSTEP=0
      Z=0.0D0
      ZE=1.0D0
      S=1.0D0
      EP=S/(1.0D4*NM)
      ZEND=ZE-EP
      DO 1 I=1,N
1     SUM(I)=(0.0D0,0.0D0)
      NS=NX
      NT=0
      CALL SAOA (Z,G1)
2     DZ=S/NS
      IF (Z+DZ.LE.ZE) GO TO 3
      DZ=ZE-Z
      IF (DZ.LE.EP) GO TO 17
3     DZOT=DZ*0.5D0
      CALL SAOA (Z+DZOT,G3)
      CALL SAOA (Z+DZ,G5)
4     NOGO=0
      DO 5 I=1,N
      T00=(G1(I)+G5(I))*DZOT
      T01(I)=(T00+DZ*G3(I))*0.5D0
      T10(I)=(4.0D0*T01(I)-T00)/3.0D0

C     Test convergence of 3 point romberg result.

      CALL TEST (DBLE(T01(I)),DBLE(T10(I)),TR,DIMAG(T01(I)),
     &           DIMAG(T10(I)),TI,0.0D0)
      IF (TR.GT.RX.OR.TI.GT.RX) NOGO=1
5     CONTINUE
      IF (NOGO.NE.0) GO TO 7
      DO 6 I=1,N
6     SUM(I)=SUM(I)+T10(I)
      NT=NT+2
      GO TO 11
7     CALL SAOA (Z+DZ*0.25D0,G2)
      CALL SAOA (Z+DZ*0.75D0,G4)
      NOGO=0
      DO 8 I=1,N
      T02=(T01(I)+DZOT*(G2(I)+G4(I)))*0.5D0
      T11=(4.0D0*T02-T01(I))/3.0D0
      T20(I)=(16.0D0*T11-T10(I))/15.0D0

C     Test convergence of 5 point romberg result.

      CALL TEST (DBLE(T11),DBLE(T20(I)),TR,DIMAG(T11),DIMAG(T20(I)),
     &           TI,0.0D0)
      IF (TR.GT.RX.OR.TI.GT.RX) NOGO=1
8     CONTINUE
      IF (NOGO.NE.0) GO TO 13
9     DO 10 I=1,N
10    SUM(I)=SUM(I)+T20(I)
      NT=NT+1
11    Z=Z+DZ
      IF (Z.GT.ZEND) GO TO 17
      DO 12 I=1,N
12    G1(I)=G5(I)
      IF (NT.LT.NTS.OR.NS.LE.NX) GO TO 2
      NS=NS/2
      NT=1
      GO TO 2
13    NT=0
      IF (NS.LT.NM) GO TO 15
      IF (LSTEP.EQ.1) GO TO 9
      LSTEP=1
      CALL LAMBDA (Z,T00,T11)
      PRINT 18, T00
      PRINT 19, Z,DZ,A,B
      DO 14 I=1,N
14    PRINT 19, G1(I),G2(I),G3(I),G4(I),G5(I)
      GO TO 9
15    NS=NS*2
      DZ=S/NS
      DZOT=DZ*0.5D0
      DO 16 I=1,N
      G5(I)=G3(I)
16    G3(I)=G2(I)
      GO TO 4
17    CONTINUE
      RETURN
C
18    FORMAT (' ROM1 -- STEP SIZE LIMITED AT LAMBDA =',2E12.5)
19    FORMAT (10E12.5)

      END
