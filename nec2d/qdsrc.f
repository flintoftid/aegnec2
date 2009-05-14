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
 
      SUBROUTINE QDSRC( IS , V , E , LD , SALP , NLOAD , NLODF , 
     &                  ZARRAY , N , M , ICON1 , ICON2 , WLAM , X , Y , 
     &                  Z , SI , BI , CAB , SAB , T1X , T1Y , T1Z , 
     &                  T2X , T2Y , T2Z , NSMAX , NQDS , IQDS , VQDS , 
     &                  JMAX , JSNO , JCO , AX , BX , CX , IFAIL , 
     &                  DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Fill incident field array for charge discontinuity voltage source. **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - IS                                                        **
C* INPUT  - V                                                         **
C* OUTPUT - E                                                         **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - NLOAD                                                     **
C* INPUT  - NLODF                                                     **
C* INPUT  - ZARRAY                                                    **
C* INPUT  - N                                                         **
C* INPUT  - M                                                         **
C* OUTPUT - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - WLAM                                                      **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - CAB                                                       **
C* INPUT  - SAB                                                       **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* INPUT  - NSMAX                                                     **
C* OUTPUT - NQDS                                                      **
C* OUTPUT - IQDS                                                      **
C* OUTPUT - VQDS                                                      **
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
C* modifies    B       CABJ    IND1    IND2    S       SABJ    SALPJ  **
C*             XJ      YJ      ZJ                                     **
C* uses value  B       CABJ    EXC     EXK     EXS     EYC     EYK    **
C*             EYS     EZC     EZK     EZS     IEXK    S       SABJ   **
C*             SALPJ                                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       EFLD    HSFLD   TBF                                    **
C* called by   ETMNS                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines
      EXTERNAL EFLD , HSFLD , TBF

C     Dummy arguments.
      INTEGER IS , LD , NLOAD , NLODF , N , M , NSMAX , NQDS , JMAX , 
     &        JSNO , IFAIL , DEBUG 
      INTEGER ICON1(2*LD) , ICON2(2*LD) , IQDS(NSMAX) , 
     &        JCO(JMAX)
      REAL*8 WLAM 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , BI(LD) , 
     &       CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX)
      COMPLEX*16 V 
      COMPLEX*16 E(1) , ZARRAY(LD) , VQDS(NSMAX)

C     Local variables.
      INTEGER I , I1 , IJ , IPR , J , JX
      REAL*8 AI , CABI , SABI , SALPI , TP , TX , TY , TZ , XI , YI , ZI
      REAL*8 CCJX(2)
      COMPLEX*16 CURD , CCJ , ETK , ETS , ETC

C     Common storage.
      INCLUDE 'dataj.inc'

C     Equivalent local storage.
      EQUIVALENCE (CCJ,CCJX)
 
C     Data initialisation.
      DATA TP /6.283185307179586476925286D0/
      DATA CCJX /0.0D0 , -0.016666666666666666666667D0/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING QDSRC'
      ENDIF
 
      I = ICON1(IS)
      ICON1(IS) = 0
      CALL TBF(IS,0,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,CX,
     &         IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      ICON1(IS) = I
      S = SI(IS)*0.5D0
      CURD = CCJ*V/((DLOG(2.0D0*S/BI(IS))-1.0D0)
     &       *(BX(JSNO)*DCOS(TP*S)+CX(JSNO)*DSIN(TP*S))*WLAM)
      NQDS = NQDS + 1
      VQDS(NQDS) = V
      IQDS(NQDS) = IS
      DO 20 JX = 1 , JSNO
         J = JCO(JX)
         S = SI(J)
         B = BI(J)
         XJ = X(J)
         YJ = Y(J)
         ZJ = Z(J)
         CABJ = CAB(J)
         SABJ = SAB(J)
         SALPJ = SALP(J)
         IF ( IEXK.EQ.0 ) GOTO 16
         IPR = ICON1(J)
         IF ( IPR ) 1 , 6 , 2
    1    IPR = -IPR
         IF ( -ICON1(IPR).NE.J ) GOTO 7
         GOTO 4
    2    IF ( IPR.NE.J ) GOTO 3
         IF ( CABJ*CABJ+SABJ*SABJ.GT.1.D-8 ) GOTO 7
         GOTO 5
    3    IF ( ICON2(IPR).NE.J ) GOTO 7
    4    XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
         IF ( XI.LT.0.999999D+0 ) GOTO 7
         IF ( DABS(BI(IPR)/B-1.0D0).GT.1.D-6 ) GOTO 7
    5    IND1 = 0
         GOTO 8
    6    IND1 = 1
         GOTO 8
    7    IND1 = 2
    8    IPR = ICON2(J)
         IF ( IPR ) 9 , 14 , 10
    9    IPR = -IPR
         IF ( -ICON2(IPR).NE.J ) GOTO 15
         GOTO 12
   10    IF ( IPR.NE.J ) GOTO 11
         IF ( CABJ*CABJ+SABJ*SABJ.GT.1.D-8 ) GOTO 15
         GOTO 13
   11    IF ( ICON1(IPR).NE.J ) GOTO 15
   12    XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
         IF ( XI.LT.0.999999D+0 ) GOTO 15
         IF ( DABS(BI(IPR)/B-1.0D0).GT.1.D-6 ) GOTO 15
   13    IND2 = 0
         GOTO 16
   14    IND2 = 1
         GOTO 16
   15    IND2 = 2
   16    CONTINUE
         DO 17 I = 1 , N
            IJ = I - J
            XI = X(I)
            YI = Y(I)
            ZI = Z(I)
            AI = BI(I)
            CALL EFLD(XI,YI,ZI,AI,IJ,IFAIL)
            IF(IFAIL.NE.0) RETURN
            CABI = CAB(I)
            SABI = SAB(I)
            SALPI = SALP(I)
            ETK = EXK*CABI + EYK*SABI + EZK*SALPI
            ETS = EXS*CABI + EYS*SABI + EZS*SALPI
            ETC = EXC*CABI + EYC*SABI + EZC*SALPI
            E(I) = E(I) - (ETK*AX(JX)+ETS*BX(JX)+ETC*CX(JX))*CURD
   17    CONTINUE
         IF ( M.EQ.0 ) GOTO 19
         IJ = LD + 1
         I1 = N
         DO 18 I = 1 , M
            IJ = IJ - 1
            XI = X(IJ)
            YI = Y(IJ)
            ZI = Z(IJ)
            CALL HSFLD(XI,YI,ZI,0.0D0)
            I1 = I1 + 1
            TX = T2X(IJ)
            TY = T2Y(IJ)
            TZ = T2Z(IJ)
            ETK = EXK*TX + EYK*TY + EZK*TZ
            ETS = EXS*TX + EYS*TY + EZS*TZ
            ETC = EXC*TX + EYC*TY + EZC*TZ
            E(I1) = E(I1) + (ETK*AX(JX)+ETS*BX(JX)+ETC*CX(JX))
     &              *CURD*SALP(IJ)
            I1 = I1 + 1
            TX = T1X(IJ)
            TY = T1Y(IJ)
            TZ = T1Z(IJ)
            ETK = EXK*TX + EYK*TY + EZK*TZ
            ETS = EXS*TX + EYS*TY + EZS*TZ
            ETC = EXC*TX + EYC*TY + EZC*TZ
            E(I1) = E(I1) + (ETK*AX(JX)+ETS*BX(JX)+ETC*CX(JX))
     &              *CURD*SALP(IJ)
   18    CONTINUE
   19    IF ( NLOAD.GT.0 .OR. NLODF.GT.0 ) E(J) = E(J) + ZARRAY(J)
     &        *CURD*(AX(JX)+CX(JX))
   20 CONTINUE
 
      RETURN
 
      END
