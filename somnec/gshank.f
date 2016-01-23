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

      SUBROUTINE GSHANK ( START , DELA , SUM , NANS , SEED , IBK , 
     &                    BK , DELB )

C*--------------------------------------------------------------------**
C*                                                                    **
C* GSHANK integrates the 6 Sommerfeld integrals from start to         **
C* infinity (until convergence) in lambda.  At the break point, BK,   **
C* the step increment may be changed from DELA to DELB.  Shank's      **
C* algorithm to accelerate convergence of a slowly converging series  **
C* is used.                                                           **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - START                                                     **
C* INPUT  - DELA                                                      **
C* OUTPUT - SUM                                                       **
C* INPUT  - NANS                                                      **
C* INPUT  - SEED                                                      **
C* INPUT  - IBK                                                       **
C* INPUT  - BK                                                        **
C* INPUT  - DELB                                                      **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    A       B                                              **
C* uses value  B                                                      **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ROM1                                                   **
C* called by   EVLUA                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL ROM1

C     Dummy arguments.
      INTEGER  NANS , IBK
      COMPLEX*16 START , DELA , SUM(6) , SEED(6) , BK , DELB

C     Local variables.
      INTEGER I , IBX , INTT , INX , J , JM , MAXH
      REAL*8 AMG , CRIT , DEN , DENM , RBK
      COMPLEX*16 A1 , A2 , AS1, AS2 , DEL , AA
      COMPLEX*16 Q1(6,20), Q2(6,20), ANS1(6), ANS2(6)

C     Save local variables.
      SAVE RBK , AMG , DEN , DENM

C     Common storage.
      INCLUDE 'cntour.inc'

C     Data initialisation.
      DATA CRIT /1.0D-4/
      DATA MAXH /20/

C      DO 99 I=1,6
C         DO 98 J=1,20
C            Q1(I,J)=(0.0D0,0.0D0)
C            Q2(I,J)=(0.0D0,0.0D0)
C 98      CONTINUE
C 99   CONTINUE
    
      RBK=DBLE(BK)
      DEL=DELA
      IBX=0
      IF (IBK.EQ.0) IBX=1
      DO 1 I=1,NANS
1     ANS2(I)=SEED(I)
      B=START
2     DO 20 INTT=1,MAXH
      INX=INTT
      A=B
      B=B+DEL
      IF (IBX.EQ.0.AND.DBLE(B).GE.RBK) GO TO 5
      CALL ROM1 (NANS,SUM,2)
      DO 3 I=1,NANS
3     ANS1(I)=ANS2(I)+SUM(I)
      A=B
      B=B+DEL
      IF (IBX.EQ.0.AND.DBLE(B).GE.RBK) GO TO 6
      CALL ROM1 (NANS,SUM,2)
      DO 4 I=1,NANS
4     ANS2(I)=ANS1(I)+SUM(I)
      GO TO 11

C     Hit break point. Reset seed and start over.

5     IBX=1
      GO TO 7
6     IBX=2
7     B=BK
      DEL=DELB
      CALL ROM1 (NANS,SUM,2)
      IF (IBX.EQ.2) GO TO 9
      DO 8 I=1,NANS
8     ANS2(I)=ANS2(I)+SUM(I)
      GO TO 2
9     DO 10 I=1,NANS
10    ANS2(I)=ANS1(I)+SUM(I)
      GO TO 2
11    DEN=0.0D0
      DO 18 I=1,NANS
      AS1=ANS1(I)
      AS2=ANS2(I)
      IF (INTT.LT.2) GO TO 17
      DO 16 J=2,INTT
      JM=J-1
      AA=Q2(I,JM)
      A1=Q1(I,JM)+AS1-2.0D0*AA
      IF (DBLE(A1).EQ.0.0D0.AND.DIMAG(A1).EQ.0.0D0) GO TO 12
      A2=AA-Q1(I,JM)
      A1=Q1(I,JM)-A2*A2/A1
      GO TO 13
12    A1=Q1(I,JM)
13    A2=AA+AS2-2.0D0*AS1
      IF (DBLE(A2).EQ.0.0D0.AND.DIMAG(A2).EQ.0.0D0) GO TO 14
      A2=AA-(AS1-AA)*(AS1-AA)/A2
      GO TO 15
14    A2=AA
15    Q1(I,JM)=AS1
      Q2(I,JM)=AS2
      AS1=A1
16    AS2=A2
17    Q1(I,INTT)=AS1
      Q2(I,INTT)=AS2
      AMG=ABS(DBLE(AS2))+DABS(DIMAG(AS2))
      IF (AMG.GT.DEN) DEN=AMG
18    CONTINUE
      DENM=1.0D-3*DEN*CRIT
      JM=INTT-3
      IF (JM.LT.1) JM=1
      DO 19 J=JM,INTT
      DO 19 I=1,NANS
      A1=Q2(I,J)
      DEN=(DABS(DBLE(A1))+DABS(DIMAG(A1)))*CRIT
      IF (DEN.LT.DENM) DEN=DENM
      A1=Q1(I,J)-A1
      AMG=ABS(DBLE(A1))+DABS(DIMAG(A1))
      IF (AMG.GT.DEN) GO TO 20
19    CONTINUE
      GO TO 22
20    CONTINUE
      PRINT 24
      DO 21 I=1,NANS
21    PRINT 25, Q1(I,INX),Q2(I,INX)
22    DO 23 I=1,NANS
23    SUM(I)=0.5D0*(Q1(I,INX)+Q2(I,INX))
      RETURN
C
24    FORMAT (' **** NO CONVERGENCE IN SUBROUTINE GSHANK ****')
25    FORMAT (10E12.5)
      END
