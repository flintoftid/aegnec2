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

      SUBROUTINE INTRP( XXX , YYY , F1 , F2 , F3 , F4 )

C*--------------------------------------------------------------------**
C*                                                                    **
C* INTRP uses bivariate cubic interpolation to obtain the values of   **
C* Sommerfeld integral contributions to the field of a source over    **
C* ground at the point (XXX,YYY).                                     **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XXX                                                       **
C* INPUT  - YYY                                                       **
C* OUTPUT - F1                                                        **
C* OUTPUT - F2                                                        **
C* OUTPUT - F3                                                        **
C* OUTPUT - F4                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ARL1    ARL2    ARL3    DXA     DYA     NXA     NYA    **
C*             XS2     XSA     YS3     YSA                            **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   SFLDS                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      REAL*8 XXX , YYY
      COMPLEX*16 F1 , F2 , F3 , F4

C     Local variables.
      INTEGER I , IADD , IADZ , IGR , IGRS , IX , IXEG , IXS , IY , 
     &        IYEG , IYS , K , ND , NDP , NXM2 , NXMS , NYM2 , NYMS
      INTEGER NDA(3) , NDPA(3)
      REAL*8 DX , DY , XXS , XS2 , XX , XZ , YYS , YS3 , YY , YZ
      COMPLEX*16 FX1 , FX2 , FX3 , FX4 , P1 , P2 , P3 , P4 , A11 , 
     &           A12 , A13 , A14 , A21 , A22 , A23 , A24 , A31 , A32 , 
     &           A33 , A34 , A41 , A42 , A43 , A44 , B11 , B12 , B13 , 
     &           B14 , B21 , B22 , B23 , B24 , B31 , B32 , B33 , B34 , 
     &           B41 , B42 , B43 , B44 , C11 , C12 , C13 , C14 , C21 , 
     &           C22 , C23 , C24 , C31 , C32 , C33 , C34 , C41 , C42 , 
     &           C43 , C44 , D11 , D12 , D13 , D14 , D21 , D22 , D23 , 
     &           D24 , D31 , D32 , D33 , D34 , D41 , D42 , D43 , D44
      COMPLEX*16 A(4,4) , BB(4,4) , C(4,4) , DD(4,4) , ARL1(1) , 
     &           ARL2(1) , ARL3(1)

C     Save local variables.
      SAVE IX , IY , IXS , IYS , IGRS , IXEG , IYEG 
      SAVE NXM2 , NYM2 , NXMS , NYMS , ND , NDP
      SAVE DX , DY , XXS , YYS , XZ , YZ     

C     Common storage.
      INCLUDE 'ggrid.inc'
 
C     Equivalent local storage.
      EQUIVALENCE (A(1,1),A11) , (A(1,2),A12) , (A(1,3),A13) , 
     &            (A(1,4),A14) , (A(2,1),A21) , (A(2,2),A22) , 
     &            (A(2,3),A23) , (A(2,4),A24) , (A(3,1),A31) , 
     &            (A(3,2),A32) , (A(3,3),A33) , (A(3,4),A34) ,
     &            (A(4,1),A41) , (A(4,2),A42) , (A(4,3),A43) , 
     &            (A(4,4),A44)
      EQUIVALENCE (BB(1,1),B11) , (BB(1,2),B12) , (BB(1,3),B13) , 
     &            (BB(1,4),B14) , (BB(2,1),B21) , (BB(2,2),B22) , 
     &            (BB(2,3),B23) , (BB(2,4),B24) , (BB(3,1),B31) , 
     &            (BB(3,2),B32) , (BB(3,3),B33) , (BB(3,4),B34) , 
     &            (BB(4,1),B41) , (BB(4,2),B42) , (BB(4,3),B43) , 
     &            (BB(4,4),B44)
      EQUIVALENCE (C(1,1),C11) , (C(1,2),C12) , (C(1,3),C13) , 
     &            (C(1,4),C14) , (C(2,1),C21) , (C(2,2),C22) , 
     &            (C(2,3),C23) , (C(2,4),C24) , (C(3,1),C31) , 
     &            (C(3,2),C32) , (C(3,3),C33) , (C(3,4),C34) ,
     &            (C(4,1),C41) , (C(4,2),C42) , (C(4,3),C43) , 
     &            (C(4,4),C44)
      EQUIVALENCE (DD(1,1),D11) , (DD(1,2),D12) , (DD(1,3),D13) , 
     &            (DD(1,4),D14) , (DD(2,1),D21) , (DD(2,2),D22) , 
     &            (DD(2,3),D23) , (DD(2,4),D24) , (DD(3,1),D31) , 
     &            (DD(3,2),D32) , (DD(3,3),D33) , (DD(3,4),D34) ,
     &            (DD(4,1),D41) , (DD(4,2),D42) , (DD(4,3),D43) , 
     &            (DD(4,4),D44)

C     Equivalent common storage.
      EQUIVALENCE (ARL1,AR1) , (ARL2,AR2) , (ARL3,AR3) , (XS2,XSA(2)) , 
     &            (YS3,YSA(3))
 
C     Data initialisation.
      DATA IXS /-10/
      DATA IYS /-10/
      DATA IGRS/-10/
      DATA DX /1.0D0/
      DATA DY /1.0D0/
      DATA XXS /0.0D0/
      DATA YYS /0.0D0/
      DATA NDA/11 , 17 , 9/
      DATA NDPA/110 , 85 , 72/
      DATA IXEG /0/
      DATA IYEG /0/
      DATA ND /0/
      DATA NDP /0/

      IF ( XXX.LT.XXS .OR. YYY.LT.YYS ) GOTO 1
      IX = INT((XXX-XXS)/DX) + 1
      IY = INT((YYY-YYS)/DY) + 1

C     If point lies in same 4 by 4 point region as previous point, old
C     values are reused.

      IF ( IX.LT.IXEG .OR. IY.LT.IYEG ) GOTO 1
      IF ( IABS(IX-IXS).LT.2 .AND. IABS(IY-IYS).LT.2 ) GOTO 13

C     Determine correct grid and grid region.

    1 IF ( XXX.GT.XS2 ) GOTO 2
      IGR = 1
      GOTO 3
    2 IGR = 2
      IF ( YYY.GT.YS3 ) IGR = 3
    3 IF ( IGR.EQ.IGRS ) GOTO 4
      IGRS = IGR
      DX = DXA(IGRS)
      DY = DYA(IGRS)
      XXS = XSA(IGRS)
      YYS = YSA(IGRS)
      NXM2 = NXA(IGRS) - 2
      NYM2 = NYA(IGRS) - 2
      NXMS = ((NXM2+1)/3)*3 + 1
      NYMS = ((NYM2+1)/3)*3 + 1
      ND = NDA(IGRS)
      NDP = NDPA(IGRS)
      IX = INT((XXX-XXS)/DX) + 1
      IY = INT((YYY-YYS)/DY) + 1
    4 IXS = ((IX-1)/3)*3 + 2
      IF ( IXS.LT.2 ) IXS = 2
      IXEG = -PATOFF
      IF ( IXS.LE.NXM2 ) GOTO 5
      IXS = NXM2
      IXEG = NXMS
    5 IYS = ((IY-1)/3)*3 + 2
      IF ( IYS.LT.2 ) IYS = 2
      IYEG = -PATOFF
      IF ( IYS.LE.NYM2 ) GOTO 6
      IYS = NYM2
      IYEG = NYMS

C     Compute coefficients of 4 cubic polynomials in XXX for the 4 grid
C     values of YYY for each of the 4 functions.

    6 IADZ = IXS + (IYS-3)*ND - NDP
      DO 12 K = 1 , 4
         IADZ = IADZ + NDP
         IADD = IADZ
         DO 11 I = 1 , 4
            IADD = IADD + ND
            GOTO (7,8,9) , IGRS
C     P1=AR1(IXS-1,IYS-2+I,K)
    7       P1 = ARL1(IADD-1)
            P2 = ARL1(IADD)
            P3 = ARL1(IADD+1)
            P4 = ARL1(IADD+2)
            GOTO 10
    8       P1 = ARL2(IADD-1)
            P2 = ARL2(IADD)
            P3 = ARL2(IADD+1)
            P4 = ARL2(IADD+2)
            GOTO 10
    9       P1 = ARL3(IADD-1)
            P2 = ARL3(IADD)
            P3 = ARL3(IADD+1)
            P4 = ARL3(IADD+2)
   10       A(I,K) = (P4-P1+3.0D0*(P2-P3))*0.1666666667D0
            BB(I,K) = (P1-2.0D0*P2+P3)*0.5D0
            C(I,K) = P3 - (2.0D0*P1+3.0D0*P2+P4)*0.1666666667D0
            DD(I,K) = P2
   11    CONTINUE
   12 CONTINUE
      XZ = (IXS-1)*DX + XXS
      YZ = (IYS-1)*DY + YYS

C     Evaluate polymomials in YYY and then use cubic interpolation in 
C     YYY for each of the 4 functions.

   13 XX = (XXX-XZ)/DX
      YY = (YYY-YZ)/DY
      FX1 = ((A11*XX+B11)*XX+C11)*XX + D11
      FX2 = ((A21*XX+B21)*XX+C21)*XX + D21
      FX3 = ((A31*XX+B31)*XX+C31)*XX + D31
      FX4 = ((A41*XX+B41)*XX+C41)*XX + D41
      P1 = FX4 - FX1 + 3.0D0*(FX2-FX3)
      P2 = 3.0D0*(FX1-2.0D0*FX2+FX3)
      P3 = 6.0D0*FX3 - 2.0D0*FX1 - 3.0D0*FX2 - FX4
      F1 = ((P1*YY+P2)*YY+P3)*YY*0.1666666667D0 + FX2
      FX1 = ((A12*XX+B12)*XX+C12)*XX + D12
      FX2 = ((A22*XX+B22)*XX+C22)*XX + D22
      FX3 = ((A32*XX+B32)*XX+C32)*XX + D32
      FX4 = ((A42*XX+B42)*XX+C42)*XX + D42
      P1 = FX4 - FX1 + 3.0D0*(FX2-FX3)
      P2 = 3.0D0*(FX1-2.0D0*FX2+FX3)
      P3 = 6.0D0*FX3 - 2.0D0*FX1 - 3.0D0*FX2 - FX4
      F2 = ((P1*YY+P2)*YY+P3)*YY*0.1666666667D0 + FX2
      FX1 = ((A13*XX+B13)*XX+C13)*XX + D13
      FX2 = ((A23*XX+B23)*XX+C23)*XX + D23
      FX3 = ((A33*XX+B33)*XX+C33)*XX + D33
      FX4 = ((A43*XX+B43)*XX+C43)*XX + D43
      P1 = FX4 - FX1 + 3.0D0*(FX2-FX3)
      P2 = 3.0D0*(FX1-2.0D0*FX2+FX3)
      P3 = 6.0D0*FX3 - 2.0D0*FX1 - 3.0D0*FX2 - FX4
      F3 = ((P1*YY+P2)*YY+P3)*YY*0.1666666667D0 + FX2
      FX1 = ((A14*XX+B14)*XX+C14)*XX + D14
      FX2 = ((A24*XX+B24)*XX+C24)*XX + D24
      FX3 = ((A34*XX+B34)*XX+C34)*XX + D34
      FX4 = ((A44*XX+B44)*XX+C44)*XX + D44
      P1 = FX4 - FX1 + 3.0D0*(FX2-FX3)
      P2 = 3.0D0*(FX1-2.0D0*FX2+FX3)
      P3 = 6.0D0*FX3 - 2.0D0*FX1 - 3.0D0*FX2 - FX4
      F4 = ((P1*YY+P2)*YY+P3)*YY*0.1666666667D0 + FX2
 
      RETURN
 
      END
