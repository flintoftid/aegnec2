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

      SUBROUTINE ROM2( A , BB , SUM , DMIN , IFAIL )

C*--------------------------------------------------------------------**
C*                                                                    **
C* For the Sommerfeld ground option, ROM2 integrates over the source  **
C* segment to obtain the total field due to ground. The method of     **
C* variable interval width Romberg integration is used. There are 9   **
C* field components - the X, Y, and ZZ components due to constant,    **
C* sine, and cosine current distributions.                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - A                                                         **
C* INPUT  - BB                                                        **
C* OUTPUT - SUM                                                       **
C* PASSED - DMIN                                                      **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       SFLDS   TEST                                           **
C* called by   EFLD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL SFLDS , TEST

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER IFAIL
      REAL*8 A , BB , DMIN
      COMPLEX*16 SUM(9)

C     Local variables.
      INTEGER I , NN , NM , NS , NT , NTS , NX
      REAL*8 DZ , DZOT , EP , RX , SSS , TI , TMAG1 , TMAG2 , TR , ZZ , 
     &       ZE , ZEND
      COMPLEX*16 T00 , T02 , T11
      COMPLEX*16 G1(9) , G2(9) , G3(9) , G4(9) , G5(9) , T01(9) ,
     &           T10(9) , T20(9)

C     Data initialisation.
      DATA NM /65536/
      DATA NTS /4/ 
      DATA NX /1/ 
      DATA NN /9/ 
      DATA RX /1.0D-4/
 
 
      ZZ = A
      ZE = BB
      SSS = BB - A
      IF ( SSS.GE.0. ) GOTO 1
      WRITE (CHRSLT,18)
      IFAIL=32
      RETURN
    1 EP = SSS/(1.0D4*NM)
      ZEND = ZE - EP
      DO 2 I = 1 , NN
         SUM(I) = (0.0D0,0.0D0)
    2 CONTINUE
      NS = NX
      NT = 0
      CALL SFLDS(ZZ,G1)
    3 DZ = SSS/NS
      IF ( ZZ+DZ.LE.ZE ) GOTO 4
      DZ = ZE - ZZ
      IF ( DZ.LE.EP ) GOTO 17
    4 DZOT = DZ*0.5D0
      CALL SFLDS(ZZ+DZOT,G3)
      CALL SFLDS(ZZ+DZ,G5)
    5 TMAG1 = 0.0D0
      TMAG2 = 0.0D0

C     Evaluate 3 point Romberg result and test convergence.

      DO 6 I = 1 , NN
         T00 = (G1(I)+G5(I))*DZOT
         T01(I) = (T00+DZ*G3(I))*0.5D0
         T10(I) = (4.0D0*T01(I)-T00)/3.0D0
         IF ( I.GT.3 ) GOTO 6
         TR = DBLE(T01(I))
         TI = DIMAG(T01(I))
         TMAG1 = TMAG1 + TR*TR + TI*TI
         TR = DBLE(T10(I))
         TI = DIMAG(T10(I))
         TMAG2 = TMAG2 + TR*TR + TI*TI
    6 CONTINUE
      TMAG1 = DSQRT(TMAG1)
      TMAG2 = DSQRT(TMAG2)
      CALL TEST(TMAG1,TMAG2,TR,0.0D0,0.0D0,TI,DMIN)
      IF ( TR.GT.RX ) GOTO 8
      DO 7 I = 1 , NN
         SUM(I) = SUM(I) + T10(I)
    7 CONTINUE
      NT = NT + 2
      GOTO 12
    8 CALL SFLDS(ZZ+DZ*0.25D0,G2)
      CALL SFLDS(ZZ+DZ*0.75D0,G4)
      TMAG1 = 0.0D0
      TMAG2 = 0.0D0

C     Evaluate 5 point Romberg result and test convergence.

      DO 9 I = 1 , NN
         T02 = (T01(I)+DZOT*(G2(I)+G4(I)))*0.5D0
         T11 = (4.0D0*T02-T01(I))/3.0D0
         T20(I) = (16.0D0*T11-T10(I))/15.0D0
         IF ( I.GT.3 ) GOTO 9
         TR = DBLE(T11)
         TI = DIMAG(T11)
         TMAG1 = TMAG1 + TR*TR + TI*TI
         TR = DBLE(T20(I))
         TI = DIMAG(T20(I))
         TMAG2 = TMAG2 + TR*TR + TI*TI
    9 CONTINUE
      TMAG1 = DSQRT(TMAG1)
      TMAG2 = DSQRT(TMAG2)
      CALL TEST(TMAG1,TMAG2,TR,0.0D0,0.0D0,TI,DMIN)
      IF ( TR.GT.RX ) GOTO 14
   10 DO 11 I = 1 , NN
         SUM(I) = SUM(I) + T20(I)
   11 CONTINUE
      NT = NT + 1
   12 ZZ = ZZ + DZ
      IF ( ZZ.GT.ZEND ) GOTO 17
      DO 13 I = 1 , NN
         G1(I) = G5(I)
   13 CONTINUE
      IF ( NT.LT.NTS .OR. NS.LE.NX ) GOTO 3
      NS = NS/2
      NT = 1
      GOTO 3
   14 NT = 0
      IF ( NS.LT.NM ) GOTO 15
      WRITE (CHRSLT,19) ZZ
      GOTO 10
   15 NS = NS*2
      DZ = SSS/NS
      DZOT = DZ*0.5D0
      DO 16 I = 1 , NN
         G5(I) = G3(I)
         G3(I) = G2(I)
   16 CONTINUE
      GOTO 5
   17 CONTINUE
 
      RETURN
 
   18 FORMAT (' ERROR - B LESS THAN A IN ROM2')
   19 FORMAT (' ROM2 -- STEP SIZE LIMITED AT Z =',1P,E12.5)
 
      END
