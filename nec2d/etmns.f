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
 
      SUBROUTINE ETMNS( P1 , P2 , P3 , P4 , P5 , P6 , IPR , E , LD , 
     &                  SALP , NLOAD , NLODF , ZARRAY , N , M , ICON1 , 
     &                  ICON2 , WLAM , X , Y , Z , SI , BI , CAB , SAB , 
     &                  T1X , T1Y , T1Z , T2X , T2Y , T2Z , NSMAX , 
     &                  NVQD , NSANT , NQDS , IVQD , ISANT , IQDS , 
     &                  VQD , VSANT , VQDS , JMAX , JSNO , JCO , AX , 
     &                  BX , CX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* ETMNS fills the array E with the negative of the electric field    **
C* incident on the structure.  E is the right hand side of the matrix **
C* equation.                                                          **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - P1                                                        **
C* INPUT  - P2                                                        **
C* INPUT  - P3                                                        **
C* INPUT  - P4                                                        **
C* INPUT  - P5                                                        **
C* INPUT  - P6                                                        **
C* INPUT  - IPR                                                       **
C* OUTPUT - E                                                         **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* PASSED - NLOAD                                                     **
C* PASSED - NLODF                                                     **
C* PASSED - ZARRAY                                                    **
C* INPUT  - N                                                         **
C* INPUT  - M                                                         **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - WLAM                                                      **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* PASSED - BI                                                        **
C* INPUT  - CAB                                                       **
C* INPUT  - SAB                                                       **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* INPUT  - NSMAX                                                     **
C* INPUT  - NVQD                                                      **
C* INPUT  - NSANT                                                     **
C* OUTPUT - NQDS                                                      **
C* INPUT  - IVQD                                                      **
C* INPUT  - ISANT                                                     **
C* PASSED - IQDS                                                      **
C* PASSED - VQD                                                       **
C* INPUT  - VSANT                                                     **
C* PASSED - VQDS                                                      **
C* INPUT  - JMAX                                                      **
C* PASSED - JSNO                                                      **
C* PASSED - JCO                                                       **
C* PASSED - AX                                                        **
C* PASSED - BX                                                        **
C* PASSED - CX                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  IPERF   KSYMP   ZRATI                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       QDSRC                                                  **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL QDSRC

C     Dummy arguments.
      INTEGER IPR , LD , NLOAD , NLODF , N , M , NSMAX , NVQD , NSANT , 
     &        NQDS , JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , IVQD(NSMAX) , 
     &        ISANT(NSMAX) , IQDS(NSMAX) , JCO(JMAX)
      REAL*8 P1 , P2 , P3 , P4 , P5 , P6 , WLAM 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , BI(LD) , 
     &       CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 E(2*LD) , ZARRAY(LD) , VQD(NSMAX) , 
     &           VSANT(NSMAX) , VQDS(NSMAX)

C     Local variables.
      INTEGER I , I1 , I2 , II , IS , NNEQ , NPM
      REAL*8 ARG , CET , CPH , CTH , DS , DSH , PX , PY , PZ , 
     &       QX , QY , QZ , R , RETA , RS , SET , SPH , STH ,
     &       TP , WX , WY , WZ , ETAOTP , ZERTOL
      COMPLEX*16 CCX , CY , CZ , ER , ET , EZH , ERH , RRV , RRH , 
     &           TT1 , TT2

C     Common storage.
      INCLUDE 'gnd.inc'
 
C     Data initialisation.
      DATA TP /6.283185307179586476925286D0/
      DATA RETA /2.654418729438072384210608D-3/
      DATA ETAOTP /59.9584916D0/
      DATA ZERTOL /1.0D-40/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING ETMNS'
      ENDIF

C     Following single line added by IDF to kill compiler warning.
      I2 = 0
C
      NNEQ = N + 2*M
      NQDS = 0
      IF ( IPR.GT.0 .AND. IPR.NE.5 ) GOTO 5

C     Applied field of voltage sources for transmitting case.

      DO 1 I = 1 , NNEQ
         E(I) = (0.0D0,0.0D0)
    1 CONTINUE
      IF ( NSANT.EQ.0 ) GOTO 3
      DO 2 I = 1 , NSANT
         IS = ISANT(I)
         E(IS) = -VSANT(I)/(SI(IS)*WLAM)
    2 CONTINUE
    3 IF ( NVQD.EQ.0 ) RETURN
      DO 4 I = 1 , NVQD
         IS = IVQD(I)
         CALL QDSRC(IS,VQD(I),E,LD,SALP,NLOAD,NLODF,ZARRAY,N,M,ICON1,
     &              ICON2,WLAM,X,Y,Z,SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,
     &              T2Z,NSMAX,NQDS,IQDS,VQDS,JMAX,JSNO,JCO,AX,BX,CX,
     &              IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
    4 CONTINUE
      RETURN
    5 IF ( IPR.GT.3 ) GOTO 19

C     Incident plane wave, linearly polarized.

      CTH = DCOS(P1)
      STH = DSIN(P1)
      CPH = DCOS(P2)
      SPH = DSIN(P2)
      CET = DCOS(P3)
      SET = DSIN(P3)
      PX = CTH*CPH*CET - SPH*SET
      PY = CTH*SPH*CET + CPH*SET
      PZ = -STH*CET
      WX = -STH*CPH
      WY = -STH*SPH
      WZ = -CTH
      QX = WY*PZ - WZ*PY
      QY = WZ*PX - WX*PZ
      QZ = WX*PY - WY*PX
      IF ( KSYMP.EQ.1 ) GOTO 7
      IF ( IPERF.EQ.1 ) GOTO 6
      RRV = CDSQRT(1.0D0-ZRATI*ZRATI*STH*STH)
      RRH = ZRATI*CTH
      RRH = (RRH-RRV)/(RRH+RRV)
      RRV = ZRATI*RRV
      RRV = -(CTH-RRV)/(CTH+RRV)
      GOTO 7
    6 RRV = -(1.0D0,0.0D0)
      RRH = -(1.0D0,0.0D0)
    7 IF ( IPR.GT.1 ) GOTO 13
      IF ( N.EQ.0 ) GOTO 10
      DO 8 I = 1 , N
         ARG = -TP*(WX*X(I)+WY*Y(I)+WZ*Z(I))
         E(I) = -(PX*CAB(I)+PY*SAB(I)+PZ*SALP(I))
     &          *DCMPLX(DCOS(ARG),DSIN(ARG))
    8 CONTINUE
      IF ( KSYMP.EQ.1 ) GOTO 10
      TT1 = (PY*CPH-PX*SPH)*(RRH-RRV)
      CCX = RRV*PX - TT1*SPH
      CY = RRV*PY + TT1*CPH
      CZ = -RRV*PZ
      DO 9 I = 1 , N
         ARG = -TP*(WX*X(I)+WY*Y(I)-WZ*Z(I))
         E(I) = E(I) - (CCX*CAB(I)+CY*SAB(I)+CZ*SALP(I))
     &          *DCMPLX(DCOS(ARG),DSIN(ARG))
    9 CONTINUE
   10 IF ( M.EQ.0 ) RETURN
      I = LD + 1
      I1 = N - 1
      DO 11 IS = 1 , M
         I = I - 1
         I1 = I1 + 2
         I2 = I1 + 1
         ARG = -TP*(WX*X(I)+WY*Y(I)+WZ*Z(I))
         TT1 = DCMPLX(DCOS(ARG),DSIN(ARG))*SALP(I)*RETA
         E(I2) = (QX*T1X(I)+QY*T1Y(I)+QZ*T1Z(I))*TT1
         E(I1) = (QX*T2X(I)+QY*T2Y(I)+QZ*T2Z(I))*TT1
   11 CONTINUE
      IF ( KSYMP.EQ.1 ) RETURN
      TT1 = (QY*CPH-QX*SPH)*(RRV-RRH)
      CCX = -(RRH*QX-TT1*SPH)
      CY = -(RRH*QY+TT1*CPH)
      CZ = RRH*QZ
      I = LD + 1
      I1 = N - 1
      DO 12 IS = 1 , M
         I = I - 1
         I1 = I1 + 2
         I2 = I1 + 1
         ARG = -TP*(WX*X(I)+WY*Y(I)-WZ*Z(I))
         TT1 = DCMPLX(DCOS(ARG),DSIN(ARG))*SALP(I)*RETA
         E(I2) = E(I2) + (CCX*T1X(I)+CY*T1Y(I)+CZ*T1Z(I))*TT1
         E(I1) = E(I1) + (CCX*T2X(I)+CY*T2Y(I)+CZ*T2Z(I))*TT1
   12 CONTINUE
      RETURN

C     Incident plane wave, elliptic polarization.

   13 TT1 = -(0.0D0,1.0D0)*P6
      IF ( IPR.EQ.3 ) TT1 = -TT1
      IF ( N.EQ.0 ) GOTO 16
      CCX = PX + TT1*QX
      CY = PY + TT1*QY
      CZ = PZ + TT1*QZ
      DO 14 I = 1 , N
         ARG = -TP*(WX*X(I)+WY*Y(I)+WZ*Z(I))
         E(I) = -(CCX*CAB(I)+CY*SAB(I)+CZ*SALP(I))
     &          *DCMPLX(DCOS(ARG),DSIN(ARG))
   14 CONTINUE
      IF ( KSYMP.EQ.1 ) GOTO 16
      TT2 = (CY*CPH-CCX*SPH)*(RRH-RRV)
      CCX = RRV*CCX - TT2*SPH
      CY = RRV*CY + TT2*CPH
      CZ = -RRV*CZ
      DO 15 I = 1 , N
         ARG = -TP*(WX*X(I)+WY*Y(I)-WZ*Z(I))
         E(I) = E(I) - (CCX*CAB(I)+CY*SAB(I)+CZ*SALP(I))
     &          *DCMPLX(DCOS(ARG),DSIN(ARG))
   15 CONTINUE
   16 IF ( M.EQ.0 ) RETURN
      CCX = QX - TT1*PX
      CY = QY - TT1*PY
      CZ = QZ - TT1*PZ
      I = LD + 1
      I1 = N - 1
      DO 17 IS = 1 , M
         I = I - 1
         I1 = I1 + 2
         I2 = I1 + 1
         ARG = -TP*(WX*X(I)+WY*Y(I)+WZ*Z(I))
         TT2 = DCMPLX(DCOS(ARG),DSIN(ARG))*SALP(I)*RETA
         E(I2) = (CCX*T1X(I)+CY*T1Y(I)+CZ*T1Z(I))*TT2
         E(I1) = (CCX*T2X(I)+CY*T2Y(I)+CZ*T2Z(I))*TT2
   17 CONTINUE
      IF ( KSYMP.EQ.1 ) RETURN
      TT1 = (CY*CPH-CCX*SPH)*(RRV-RRH)
      CCX = -(RRH*CCX-TT1*SPH)
      CY = -(RRH*CY+TT1*CPH)
      CZ = RRH*CZ
      I = LD + 1
      I1 = N - 1
      DO 18 IS = 1 , M
         I = I - 1
         I1 = I1 + 2
         I2 = I1 + 1
         ARG = -TP*(WX*X(I)+WY*Y(I)-WZ*Z(I))
         TT1 = DCMPLX(DCOS(ARG),DSIN(ARG))*SALP(I)*RETA
         E(I2) = E(I2) + (CCX*T1X(I)+CY*T1Y(I)+CZ*T1Z(I))*TT1
         E(I1) = E(I1) + (CCX*T2X(I)+CY*T2Y(I)+CZ*T2Z(I))*TT1
   18 CONTINUE
      RETURN

C     Incident field of an elementary current source.

   19 WZ = DCOS(P4)
      WX = WZ*DCOS(P5)
      WY = WZ*DSIN(P5)
      WZ = DSIN(P4)
      DS = P6*ETAOTP
      DSH = P6/(2.0D0*TP)
      NPM = N + M
      IS = LD + 1
      I1 = N - 1
      DO 24 I = 1 , NPM
         II = I
         IF ( I.LE.N ) GOTO 20
         IS = IS - 1
         II = IS
         I1 = I1 + 2
         I2 = I1 + 1
   20    PX = X(II) - P1
         PY = Y(II) - P2
         PZ = Z(II) - P3
         RS = PX*PX + PY*PY + PZ*PZ
         IF ( RS.LT.ZERTOL ) GOTO 24
         R = DSQRT(RS)
         PX = PX/R
         PY = PY/R
         PZ = PZ/R
         CTH = PX*WX + PY*WY + PZ*WZ
         STH = DSQRT(1.0D0-CTH*CTH)
         QX = PX - WX*CTH
         QY = PY - WY*CTH
         QZ = PZ - WZ*CTH
         ARG = DSQRT(QX*QX+QY*QY+QZ*QZ)
         IF ( ARG.LT.ZERTOL ) GOTO 21
         QX = QX/ARG
         QY = QY/ARG
         QZ = QZ/ARG
         GOTO 22
   21    QX = 1.0D0
         QY = 0.0D0
         QZ = 0.0D0
   22    ARG = -TP*R
         TT1 = DCMPLX(DCOS(ARG),DSIN(ARG))
         IF ( I.GT.N ) GOTO 23
         TT2 = DCMPLX(1.0D0,-1.0D0/(R*TP))/RS
         ER = DS*TT1*TT2*CTH
         ET = 0.5D0*DS*TT1*((0.0D0,1.0D0)*TP/R+TT2)*STH
         EZH = ER*CTH - ET*STH
         ERH = ER*STH + ET*CTH
         CCX = EZH*WX + ERH*QX
         CY = EZH*WY + ERH*QY
         CZ = EZH*WZ + ERH*QZ
         E(I) = -(CCX*CAB(I)+CY*SAB(I)+CZ*SALP(I))
         GOTO 24
   23    PX = WY*QZ - WZ*QY
         PY = WZ*QX - WX*QZ
         PZ = WX*QY - WY*QX
         TT2 = DSH*TT1*DCMPLX(1.0D0/R,TP)/R*STH*SALP(II)
         CCX = TT2*PX
         CY = TT2*PY
         CZ = TT2*PZ
         E(I2) = CCX*T1X(II) + CY*T1Y(II) + CZ*T1Z(II)
         E(I1) = CCX*T2X(II) + CY*T2Y(II) + CZ*T2Z(II)
   24 CONTINUE
 
      RETURN
 
      END
