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

      SUBROUTINE EVLUA ( ERV , EZV , ERH , EPH )

C*--------------------------------------------------------------------**
C*                                                                    **
C* EVLUA controls the integration contour in the complex lambda       **
C* plane for evaluation of the Sommerfeld integrals.                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - ERV                                                       **
C* OUTPUT - EZV                                                       **
C* OUTPUT - ERH                                                       **
C* OUTPUT - EPH                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    A       B       JH                                     **
C* uses value  B       CK1     CK1SQ   CK2     CK2SQ   RHO     TKMAG  **
C*             ZPH                                                    **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       GSHANK  ROM1                                           **
C* called by   SOMNEC                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL GSHANK , ROM1

C     Dummy arguments.
      COMPLEX*16 ERV , EZV , ERH , EPH

C     Local variables.
      INTEGER I
      REAL*8 DEL , PTP , RMIS , SLOPE
      COMPLEX*16 BK , DELTA , DELTA2 , CP1 , CP2 , CP3
      COMPLEX*16 SUM(6), ANS(6)

C     Save local variables.
      SAVE DEL , SLOPE , RMIS , CP1 , CP2 , CP3 , BK , DELTA , DELTA2
C     SAVE SUM , ANS

C     Common storage.
      INCLUDE 'cntour.inc'
      INCLUDE 'evlcom.inc'

C     Data initialisation.
      DATA PTP /0.6283185308D0/

      BK=(0.0D0,0.0D0)

      DEL=ZPH
      IF (RHO.GT.DEL) DEL=RHO
      IF (ZPH.LT.2.0D0*RHO) GO TO 4

C     Bessel function form of Sommerfeld integrals.

      JH=0
      A=(0.0D0,0.0D0)
      DEL=1.0D0/DEL
      IF (DEL.LE.TKMAG) GO TO 2
      B=DCMPLX(0.1D0*TKMAG,-0.1D0*TKMAG)
      CALL ROM1 (6,SUM,2)
      A=B
      B=DCMPLX(DEL,-DEL)
      CALL ROM1 (6,ANS,2)
      DO 1 I=1,6
1     SUM(I)=SUM(I)+ANS(I)
      GO TO 3
2     B=DCMPLX(DEL,-DEL)
      CALL ROM1 (6,SUM,2)
3     DELTA=PTP*DEL
      CALL GSHANK (B,DELTA,ANS,6,SUM,0,B,B)
      GO TO 10

C     Hankel function form of Sommerfeld integrals.

4     JH=1
      CP1=DCMPLX(0.0D0,0.4D0*CK2)
      CP2=DCMPLX(0.6D0*CK2,-0.2D0*CK2)
      CP3=DCMPLX(1.02D0*CK2,-0.2D0*CK2)
      A=CP1
      B=CP2
      CALL ROM1 (6,SUM,2)
      A=CP2
      B=CP3
      CALL ROM1 (6,ANS,2)
      DO 5 I=1,6
5     SUM(I)=-(SUM(I)+ANS(I))

C     Path from imaginary axis to -infinity.

      SLOPE=1000.0D0
      IF (ZPH.GT.0.001D0*RHO) SLOPE=RHO/ZPH
      DEL=PTP/DEL
      DELTA=DCMPLX(-1.0D0,SLOPE)*DEL/SQRT(1.0D0+SLOPE*SLOPE)
      DELTA2=-DCONJG(DELTA)
      CALL GSHANK (CP1,DELTA,ANS,6,SUM,0,BK,BK)
      RMIS=RHO*(DBLE(CK1)-CK2)
      IF (RMIS.LT.2.0D0*CK2) GO TO 8
      IF (RHO.LT.1.0D-10) GO TO 8
      IF (ZPH.LT.1.0D-10) GO TO 6
      BK=DCMPLX(-ZPH,RHO)*(CK1-CP3)
      RMIS=-DBLE(BK)/ABS(DIMAG(BK))
      IF(RMIS.GT.4.0D0*RHO/ZPH)GO TO 8

C     Integrate up between branch cuts, then to + infinity.

6     CP1=CK1-(0.1D0,0.2D0)
      CP2=CP1+0.2D0
      BK=DCMPLX(0.0D0,DEL)
      CALL GSHANK (CP1,BK,SUM,6,ANS,0,BK,BK)
      A=CP1
      B=CP2
      CALL ROM1 (6,ANS,1)
      DO 7 I=1,6
7     ANS(I)=ANS(I)-SUM(I)
      CALL GSHANK (CP3,BK,SUM,6,ANS,0,BK,BK)
      CALL GSHANK (CP2,DELTA2,ANS,6,SUM,0,BK,BK)
      GO TO 10

C     Integrate below branch points, then to + infinity.

8     DO 9 I=1,6
9     SUM(I)=-ANS(I)
      RMIS=DBLE(CK1)*1.01D0
      IF (CK2+1.0D0.GT.RMIS) RMIS=CK2+1.0D0
      BK=DCMPLX(RMIS,0.99D0*DIMAG(CK1))
      DELTA=BK-CP3
      DELTA=DELTA*DEL/CDABS(DELTA)
      CALL GSHANK (CP3,DELTA,ANS,6,SUM,1,BK,DELTA2)
10    ANS(6)=ANS(6)*CK1

C     Conjugate since nec uses exp(+jwt).

      ERV=DCONJG(CK1SQ*ANS(3))
      EZV=DCONJG(CK1SQ*(ANS(2)+CK2SQ*ANS(5)))
      ERH=DCONJG(CK2SQ*(ANS(1)+ANS(6)))
      EPH=-DCONJG(CK2SQ*(ANS(4)+ANS(6)))

      RETURN

      END
