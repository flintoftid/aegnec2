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
 
      SUBROUTINE CABC( CURX , LD , AIR , AII , BIR , BII , CIR , 
     &                 CII , N , M , ICON1 , ICON2 , WLAM , SI , BI , 
     &                 T1X , T1Y , T1Z , T2X , T2Y , T2Z , NSMAX , 
     &                 NQDS , IQDS , VQDS , JMAX , JSNO , JCO , AX , 
     &                 BX , CX , IFAIL , DEBUG )
 
C*--------------------------------------------------------------------**
C*                                                                    **
C* CABC computes coefficients of the constant (A), sine (B), and      **
C* cosine (C) terms in the current interpolation functions for the    **
C* current vector CUR. Surface current components are also computed.  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - CURX                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - AIR                                                       **
C* OUTPUT - AII                                                       **
C* OUTPUT - BIR                                                       **
C* OUTPUT - BII                                                       **
C* OUTPUT - CIR                                                       **
C* OUTPUT - CII                                                       **
C* INPUT  - N                                                         **
C* INPUT  - M                                                         **
C* OUTPUT - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - WLAM                                                      **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* INPUT  - NSMAX                                                     **
C* INPUT  - NQDS                                                      **
C* INPUT  - IQDS                                                      **
C* INPUT  - VQDS                                                      **
C* INPUT  - JMAX                                                      **
C* INPUT  - JSNO                                                      **
C* INPUT  - JCO                                                       **
C* INPUT  - AX                                                        **
C* INPUT  - BX                                                        **
C* INPUT  - CX                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       TBF                                                    **
C* called by   NETWK                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL TBF

C     Dummy arguments.
      INTEGER LD , N , M , NSMAX , NQDS , JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , IQDS(NSMAX) , 
     &        JCO(JMAX)
      REAL*8 WLAM 
      REAL*8 AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , SI(LD) , BI(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 CURX(1) , VQDS(NSMAX)

C     Local variables.
      INTEGER I , IS , J , JCO1 , JCO2 , JX , K , JXX
      REAL*8 AI , AR , SH , TP , CCJX(2)
      COMPLEX*16 CURD , CCJ , CS1 , CS2

C     Equivalent local storage.
      EQUIVALENCE (CCJ,CCJX)

C     Data initialisation. 
      DATA TP/6.283185307179586476925286D0/
      DATA CCJX/0.0D0 , -0.016666666666666666666667D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEDUG: ENTERING CABC'
      ENDIF

      IF ( N.EQ.0 ) GOTO 8
      DO 1 I = 1 , N
         AIR(I) = 0.0D0
         AII(I) = 0.0D0
         BIR(I) = 0.0D0
         BII(I) = 0.0D0
         CIR(I) = 0.0D0
         CII(I) = 0.0D0
    1 CONTINUE
      DO 3 I = 1 , N
         AR = DBLE(CURX(I))
         AI = DIMAG(CURX(I))
         CALL TBF(I,1,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,CX,
     &            IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         DO 2 JX = 1 , JSNO
            J = JCO(JX)
            AIR(J) = AIR(J) + AX(JX)*AR
            AII(J) = AII(J) + AX(JX)*AI
            BIR(J) = BIR(J) + BX(JX)*AR
            BII(J) = BII(J) + BX(JX)*AI
            CIR(J) = CIR(J) + CX(JX)*AR
            CII(J) = CII(J) + CX(JX)*AI
    2    CONTINUE
    3 CONTINUE
      IF ( NQDS.EQ.0 ) GOTO 6
      DO 5 IS = 1 , NQDS
         I = IQDS(IS)
         JXX = ICON1(I)
         ICON1(I) = 0
         CALL TBF(I,0,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,CX,
     &            IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         ICON1(I) = JXX
         SH = SI(I)*0.5D0
         CURD = CCJ*VQDS(IS)
     &          /((DLOG(2.0D0*SH/BI(I))-1.0D0)
     &          *(BX(JSNO)*DCOS(TP*SH)+CX(JSNO)
     &          *DSIN(TP*SH))*WLAM)
         AR = DBLE(CURD)
         AI = DIMAG(CURD)
         DO 4 JX = 1 , JSNO
            J = JCO(JX)
            AIR(J) = AIR(J) + AX(JX)*AR
            AII(J) = AII(J) + AX(JX)*AI
            BIR(J) = BIR(J) + BX(JX)*AR
            BII(J) = BII(J) + BX(JX)*AI
            CIR(J) = CIR(J) + CX(JX)*AR
            CII(J) = CII(J) + CX(JX)*AI
    4    CONTINUE
    5 CONTINUE
    6 DO 7 I = 1 , N
         CURX(I) = DCMPLX(AIR(I)+CIR(I),AII(I)+CII(I))
    7 CONTINUE
    8 IF ( M.EQ.0 ) RETURN
C     Convert surface currents from T1,t2 components to X,Y,Z 
C     components.
      K = LD - M
      JCO1 = N + 2*M + 1
      JCO2 = JCO1 + M
      DO 9 I = 1 , M
         K = K + 1
         JCO1 = JCO1 - 2
         JCO2 = JCO2 - 3
         CS1 = CURX(JCO1)
         CS2 = CURX(JCO1+1)
         CURX(JCO2) = CS1*T1X(K) + CS2*T2X(K)
         CURX(JCO2+1) = CS1*T1Y(K) + CS2*T2Y(K)
         CURX(JCO2+2) = CS1*T1Z(K) + CS2*T2Z(K)
    9 CONTINUE
 
      RETURN
 
      END
