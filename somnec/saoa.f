C* 
C* aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
C* Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>
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

      SUBROUTINE SAOA ( T , ANS )

C*--------------------------------------------------------------------**
C*                                                                    **
C* SAOA computes the integrand for each of the 6 Sommerfeld integrals **
C* for source and observer above ground.                              **
C*                                                                    **
**------------------------- Dummy Arguments --------------------------**
**                                                                    **
*A PASSED - T                                                         **
*A OUTPUT - ANS                                                       **
**                                                                    **
**------------------------- COMMON Variables -------------------------**
**                                                                    **
** modifies    ** NOTHING **                                          **
** uses value  CK1     CK1R    CK1SQ   CK2     CK2SQ   CKSM    CT1    **
**             CT2     CT3     JH      RHO     TSMAG   ZPH            **
**                                                                    **
**----------------------- External Subprograms -----------------------**
**                                                                    **
** calls       BESSEL  HANKEL  LAMBDA                                 **
** called by   ROM1                                                   **
**                                                                    **
**--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL BESSEL , HANKEL , LAMBDA

C     Dummy arguments.
      REAL*8 T
      COMPLEX*16 ANS(6)

C     Local variables.
      REAL*8 SIGN , XLR
      COMPLEX*16 XL , DXL , CGAM1 , CGAM2 , B0 , B0P , COM , DGAM , 
     &           DEN1 , DEN2

C     Save local variables.
      SAVE XL , DXL , CGAM1 , CGAM2 , B0 , B0P , COM , DGAM , DEN1 , 
     &     DEN2

C     Common storage.
      INCLUDE 'evlcom.inc'

      CALL LAMBDA (T,XL,DXL)
      IF (JH.GT.0) GO TO 1

C     Bessel function form.

      CALL BESSEL (XL*RHO,B0,B0P)
      B0=2.0D0*B0
      B0P=2.0D0*B0P
      CGAM1=CDSQRT(XL*XL-CK1SQ)
      CGAM2=CDSQRT(XL*XL-CK2SQ)
      IF (DBLE(CGAM1).EQ.0.0D0) CGAM1=DCMPLX(0.0D0,-DABS(DIMAG(CGAM1)))
      IF (DBLE(CGAM2).EQ.0.0D0) CGAM2=DCMPLX(0.0D0,-DABS(DIMAG(CGAM2)))
      GO TO 2

C     Hankel function form.

1     CALL HANKEL (XL*RHO,B0,B0P)
      COM=XL-CK1
      CGAM1=CDSQRT(XL+CK1)*CDSQRT(COM)
      IF (DBLE(COM).LT.0.0D0.AND.DIMAG(COM).GE.0.0D0) CGAM1=-CGAM1
      COM=XL-CK2
      CGAM2=CDSQRT(XL+CK2)*CDSQRT(COM)
      IF (DBLE(COM).LT.0.0D0.AND.DIMAG(COM).GE.0.0D0) CGAM2=-CGAM2
2     XLR=DBLE(XL*DCONJG(XL))
      IF (XLR.LT.TSMAG) GO TO 3
      IF (DIMAG(XL).LT.0.0D0) GO TO 4
      XLR=DBLE(XL)
      IF (XLR.LT.CK2) GO TO 5
      IF (XLR.GT.CK1R) GO TO 4
3     DGAM=CGAM2-CGAM1
      GO TO 7
4     SIGN=1.0D0
      GO TO 6
5     SIGN=-1.0D0
6     DGAM=1.0D0/(XL*XL)
      DGAM=SIGN*((CT3*DGAM+CT2)*DGAM+CT1)/XL
7     DEN2=CKSM*DGAM/(CGAM2*(CK1SQ*CGAM2+CK2SQ*CGAM1))
      DEN1=1.0D0/(CGAM1+CGAM2)-CKSM/CGAM2
      COM=DXL*XL*CDEXP(-CGAM2*ZPH)
      ANS(6)=COM*B0*DEN1/CK1
      COM=COM*DEN2
      IF (RHO.EQ.0.0D0) GO TO 8
      B0P=B0P/RHO
      ANS(1)=-COM*XL*(B0P+B0*XL)
      ANS(4)=COM*XL*B0P
      GO TO 9
8     ANS(1)=-COM*XL*XL*0.5D0
      ANS(4)=ANS(1)
9     ANS(2)=COM*CGAM2*CGAM2*B0
      ANS(3)=-ANS(4)*CGAM2*RHO
      ANS(5)=COM*B0
      RETURN

      END
