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

      SUBROUTINE CMWS( J , I1 , I2 , CCM , NR , CW , NW , ITRP , 
     &                 LD , SALP , X , Y , Z , SI , BI , CAB , 
     &                 SAB , T1X , T1Y , T1Z , T2X , T2Y , T2Z , JMAX , 
     &                 JSNO , JCO , AX , BX , CX , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CMWS computes and stores matrix elements for the H field at patch  **
C* centres due to current on wire segments.                           **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - J                                                         **
C* INPUT  - I1                                                        **
C* INPUT  - I2                                                        **
C* OUTPUT - CCM                                                       **
C* INPUT  - NR                                                        **
C* OUTPUT - CW                                                        **
C* INPUT  - NW                                                        **
C* INPUT  - ITRP                                                      **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
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
C* INPUT  - JMAX                                                      **
C* INPUT  - JSNO                                                      **
C* INPUT  - JCO                                                       **
C* INPUT  - AX                                                        **
C* INPUT  - BX                                                        **
C* INPUT  - CX                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    B       CABJ    S       SABJ    SALPJ   XJ      YJ     **
C*             ZJ                                                     **
C* uses value  EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS                                            **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       HSFLD                                                  **
C* called by   CMNGF   CMSET                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL HSFLD

C     Dummy arguments.
      INTEGER J , I1 , I2 , NR , NW , ITRP , LD , JMAX , JSNO , DEBUG
      INTEGER JCO(JMAX)
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , CAB(LD) , SAB(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 CCM(NR,1) , CW(NW,1)

C     Local variables.
      INTEGER I , IJ , IK , IPATCH , IPR , JS , JX , LDP
      REAL*8 TX , TY , TZ , XI , YI , ZI
      COMPLEX*16 ETK , ETS , ETC

C     Common storage.
      INCLUDE 'dataj.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMWS'
      ENDIF
 
C     Following one line added by IDF to kill compiler warning.
      JS = 0
C
      LDP = LD + 1
      S = SI(J)
      B = BI(J)
      XJ = X(J)
      YJ = Y(J)
      ZJ = Z(J)
      CABJ = CAB(J)
      SABJ = SAB(J)
      SALPJ = SALP(J)

C     Observation loop.

      IPR = 0
      DO 9 I = I1 , I2
         IPR = IPR + 1
         IPATCH = (I+1)/2
         IK = I - (I/2)*2
         IF ( IK.EQ.0 .AND. IPR.NE.1 ) GOTO 1
         JS = LDP - IPATCH
         XI = X(JS)
         YI = Y(JS)
         ZI = Z(JS)
         CALL HSFLD(XI,YI,ZI,0.0D0)
         IF ( IK.EQ.0 ) GOTO 1
         TX = T2X(JS)
         TY = T2Y(JS)
         TZ = T2Z(JS)
         GOTO 2
    1    TX = T1X(JS)
         TY = T1Y(JS)
         TZ = T1Z(JS)
    2    ETK = -(EXK*TX+EYK*TY+EZK*TZ)*SALP(JS)
         ETS = -(EXS*TX+EYS*TY+EZS*TZ)*SALP(JS)
         ETC = -(EXC*TX+EYC*TY+EZC*TZ)*SALP(JS)

C     Fill matrix elements. Element locations determined by connection
C     data.

         IF ( ITRP.NE.0 ) GOTO 4
C     Normal fill.
         DO 3 IJ = 1 , JSNO
            JX = JCO(IJ)
            CCM(IPR,JX) = CCM(IPR,JX) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
    3    CONTINUE
         GOTO 9
    4    IF ( ITRP.EQ.2 ) GOTO 6
C     Transposed fill.
         DO 5 IJ = 1 , JSNO
            JX = JCO(IJ)
            CCM(JX,IPR) = CCM(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
    5    CONTINUE
         GOTO 9
C     Transposed fill - C(WS) and D(WS)prime (=CW).
    6    DO 8 IJ = 1 , JSNO
            JX = JCO(IJ)
            IF ( JX.GT.NR ) GOTO 7
            CCM(JX,IPR) = CCM(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
            GOTO 8
    7       JX = JX - NR
            CW(JX,IPR) = CW(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
    8    CONTINUE
    9 CONTINUE
 
      RETURN
 
      END
