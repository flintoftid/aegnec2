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
**--------------------------------------------------------------------**
 
      SUBROUTINE CMWW( J , I1 , I2 , CCM , NR , CW , NW , ITRP , 
     &                 LD , SALP , ICON1 , ICON2 , X , Y , Z , SI , 
     &                 BI , CAB , SAB , JMAX , JSNO , JCO , AX , BX , 
     &                 CX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CMWW computes and stores matrix elements for the E field at        **
C* segment centres due to current on other segments.                  **
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
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - CAB                                                       **
C* INPUT  - SAB                                                       **
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
C*             EYS     EZC     EZK     EZS     IEXK    SABJ    SALPJ  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       EFLD                                                   **
C* called by   CMNGF   CMSET                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL EFLD

C     Dummy arguments.
      INTEGER J , I1 , I2 , NR , NW , ITRP , LD , JMAX , JSNO , IFAIL ,
     &        DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , JCO(JMAX) 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , CAB(LD) , SAB(LD) , 
     &       AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 CCM(NR,1) , CW(NW,1)

C     Local variables.
      INTEGER I , IJ , IPR , JX 
      REAL*8 AI , CABI , SABI , SALPI , XI , YI , ZI, COLTST
      COMPLEX*16 ETK , ETS , ETC

C     Common storage.
      INCLUDE 'dataj.inc'

C     Data initialisation.
      DATA COLTST /0.999999D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMWW'
      ENDIF
 
C     Set source segment parameters.
      S = SI(J)
      B = BI(J)
      XJ = X(J)
      YJ = Y(J)
      ZJ = Z(J)
      CABJ = CAB(J)
      SABJ = SAB(J)
      SALPJ = SALP(J)
      IF ( IEXK.EQ.0 ) GOTO 16
C     Decide wether ext. T.W. approx. can be used.
      IPR = ICON1(J)
      IF ( IPR ) 1 , 6 , 2
    1 IPR = -IPR
      IF ( -ICON1(IPR).NE.J ) GOTO 7
      GOTO 4
    2 IF ( IPR.NE.J ) GOTO 3
      IF ( CABJ*CABJ+SABJ*SABJ.GT.1.0D-8 ) GOTO 7
      GOTO 5
    3 IF ( ICON2(IPR).NE.J ) GOTO 7
    4 XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
      IF ( XI.LT.COLTST ) GOTO 7
      IF ( DABS(BI(IPR)/B-1.0D0).GT.1.0D-6 ) GOTO 7
    5 IND1 = 0
      GOTO 8
    6 IND1 = 1
      GOTO 8
    7 IND1 = 2
    8 IPR = ICON2(J)
      IF ( IPR ) 9 , 14 , 10
    9 IPR = -IPR
      IF ( -ICON2(IPR).NE.J ) GOTO 15
      GOTO 12
   10 IF ( IPR.NE.J ) GOTO 11
      IF ( CABJ*CABJ+SABJ*SABJ.GT.1.0D-8 ) GOTO 15
      GOTO 13
   11 IF ( ICON1(IPR).NE.J ) GOTO 15
   12 XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
      IF ( XI.LT.COLTST ) GOTO 15
      IF ( DABS(BI(IPR)/B-1.0D0).GT.1.0D-6 ) GOTO 15
   13 IND2 = 0
      GOTO 16
   14 IND2 = 1
      GOTO 16
   15 IND2 = 2
   16 CONTINUE

C     Observation loop.

      IPR = 0
      DO 23 I = I1 , I2
         IPR = IPR + 1
         IJ = I - J
         XI = X(I)
         YI = Y(I)
         ZI = Z(I)
         AI = BI(I)
         CABI = CAB(I)
         SABI = SAB(I)
         SALPI = SALP(I)
         CALL EFLD(XI,YI,ZI,AI,IJ,IFAIL)
         IF(IFAIL.NE.0) RETURN
         ETK = EXK*CABI + EYK*SABI + EZK*SALPI
         ETS = EXS*CABI + EYS*SABI + EZS*SALPI
         ETC = EXC*CABI + EYC*SABI + EZC*SALPI

C     Fill matrix elements. Element locations determined by connection
C     data.

         IF ( ITRP.NE.0 ) GOTO 18
C     Normal fill.
         DO 17 IJ = 1 , JSNO
            JX = JCO(IJ)
            CCM(IPR,JX) = CCM(IPR,JX) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
   17    CONTINUE
         GOTO 23
   18    IF ( ITRP.EQ.2 ) GOTO 20
C     Transposed fill.
         DO 19 IJ = 1 , JSNO
            JX = JCO(IJ)
            CCM(JX,IPR) = CCM(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
   19    CONTINUE
         GOTO 23
C     Trans. fill for C(WW) - test for elements for D(WW)prime (=CW).
   20    DO 22 IJ = 1 , JSNO
            JX = JCO(IJ)
            IF ( JX.GT.NR ) GOTO 21
            CCM(JX,IPR) = CCM(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
            GOTO 22
   21       JX = JX - NR
            CW(JX,IPR) = CW(JX,IPR) + ETK*AX(IJ) + ETS*BX(IJ)
     &                   + ETC*CX(IJ)
   22    CONTINUE
   23 CONTINUE
 
      RETURN
 
      END
