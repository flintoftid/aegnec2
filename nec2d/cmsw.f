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
**--------------------------------------------------------------------**

      SUBROUTINE CMSW( J1 , J2 , I1 , I2 , CCM , CW , NCW , NROW , 
     &                 ITRP , LD , SALP , N1 , N ,NP , M1 , M , MP , 
     &                 ICON1 , ICON2 , ICONX , X , Y , Z , SI , BI , 
     &                 CAB , SAB , T1X , T1Y , T1Z , T2X , T2Y , T2Z , 
     &                 JMAX , JSNO , JCO , AX , BX , CX , IFAIL , 
     &                 DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Computes and stores matrix elements for the E field at segment     **
C* centres due to current on patches.                                 **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - J1                                                        **
C* INPUT  - J2                                                        **
C* INPUT  - I1                                                        **
C* INPUT  - I2                                                        **
C* OUTPUT - CCM                                                       **
C* OUTPUT - CW                                                        **
C* INPUT  - NCW                                                       **
C* INPUT  - NROW                                                      **
C* INPUT  - ITRP                                                      **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - N1                                                        **
C* INPUT  - N                                                         **
C* INPUT  - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M                                                         **
C* INPUT  - MP                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - ICONX                                                     **
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
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXC     IPGND   S       T1XJ    T1YJ    T1ZJ    T2XJ   **
C*             T2YJ    T2ZJ    XJ      YJ      ZJ                     **
C* uses value  EXC     EXK     EXS     EYK     EYS     EZK     EZS    **
C*             KSYMP                                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       PCINT   TRIO    UNERE                                  **
C* called by   CMNGF   CMSET                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL PCINT , TRIO , UNERE

C     Parameter defintions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER J1 , J2 , I1 , I2 , NCW , NROW , ITRP , LD , N1 , N , NP , 
     &        M1 , M , MP , JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , ICONX(LD) , 
     &        JCO(JMAX)
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , CAB(LD) , SAB(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) ,AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 CCM(NROW,1) , CW(NROW,1)

C     Local variables.
      INTEGER I , ICGO , IL , IIP , IPCH , J , JL , JS , K , LDP , NEQS
      REAL*8 CABI , FSIGN , PI , PX , PY , SABI , SALPI , XI , YI , ZI
      COMPLEX*16 EMEL(9)

C     Common storage.
      INCLUDE 'gnd.inc'
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
 
C     Data initialisation.
      DATA PI/3.141592653589793238462643D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMSW'
      ENDIF

C     Following one line addedc by idf in response to compile warning.
      FSIGN=0.0D0
C
      LDP = LD + 1
      NEQS = N - N1 + 2*(M-M1)
      IF ( ITRP.LT.0 ) GOTO 15
      K = 0
      ICGO = 1
C     Observation loop.
      DO 14 I = I1 , I2
         K = K + 1
         XI = X(I)
         YI = Y(I)
         ZI = Z(I)
         CABI = CAB(I)
         SABI = SAB(I)
         SALPI = SALP(I)
         IPCH = 0
         IF ( ICON1(I).LT.PATOFF ) GOTO 1
         IPCH = ICON1(I) - PATOFF
         FSIGN = -1.0D0
    1    IF ( ICON2(I).LT.PATOFF ) GOTO 2
         IPCH = ICON2(I) - PATOFF
         FSIGN = 1.0D0
    2    JL = 0
C     Source loop.
         DO 13 J = J1 , J2
            JS = LDP - J
            JL = JL + 2
            T1XJ = T1X(JS)
            T1YJ = T1Y(JS)
            T1ZJ = T1Z(JS)
            T2XJ = T2X(JS)
            T2YJ = T2Y(JS)
            T2ZJ = T2Z(JS)
            XJ = X(JS)
            YJ = Y(JS)
            ZJ = Z(JS)
            S = BI(JS)
C     Ground loop.
            DO 12 IIP = 1 , KSYMP
               IPGND = IIP
               IF ( IPCH.NE.J .AND. ICGO.EQ.1 ) GOTO 9
               IF ( IIP.EQ.2 ) GOTO 9
               IF ( ICGO.GT.1 ) GOTO 6
               CALL PCINT(XI,YI,ZI,CABI,SABI,SALPI,EMEL)
               PY = PI*SI(I)*FSIGN
               PX = DSIN(PY)
               PY = DCOS(PY)
               EXC = EMEL(9)*FSIGN
               CALL TRIO(I,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                   CX,IFAIL,DEBUG)
               IF(IFAIL.EQ.0) RETURN
               IF ( I.GT.N1 ) GOTO 3
               IL = NEQS + ICONX(I)
               GOTO 4
    3          IL = I - NCW
               IF ( I.LE.NP ) IL = ((IL-1)/NP)*2*MP + IL
    4          IF ( ITRP.NE.0 ) GOTO 5
               CW(K,IL) = CW(K,IL)
     &                    + EXC*(AX(JSNO)+BX(JSNO)*PX+CX(JSNO)*PY)
               GOTO 6
    5          CW(IL,K) = CW(IL,K)
     &                    + EXC*(AX(JSNO)+BX(JSNO)*PX+CX(JSNO)*PY)
    6          IF ( ITRP.NE.0 ) GOTO 7
               CCM(K,JL-1) = EMEL(ICGO)
               CCM(K,JL) = EMEL(ICGO+4)
               GOTO 8
    7          CCM(JL-1,K) = EMEL(ICGO)
               CCM(JL,K) = EMEL(ICGO+4)
    8          ICGO = ICGO + 1
               IF ( ICGO.EQ.5 ) ICGO = 1
               GOTO 11
    9          CALL UNERE(XI,YI,ZI)
               IF ( ITRP.NE.0 ) GOTO 10
C     Normal fill.
               CCM(K,JL-1) = CCM(K,JL-1) + EXK*CABI + EYK*SABI 
     &                     + EZK*SALPI
               CCM(K,JL) = CCM(K,JL) + EXS*CABI + EYS*SABI 
     &                   + EZS*SALPI
               GOTO 11
C     Transposed fill.
   10          CCM(JL-1,K) = CCM(JL-1,K) + EXK*CABI + EYK*SABI 
     &                     + EZK*SALPI
               CCM(JL,K) = CCM(JL,K) + EXS*CABI + EYS*SABI 
     &                   + EZS*SALPI
   11          CONTINUE
   12       CONTINUE
   13    CONTINUE
   14 CONTINUE
      RETURN
C     For old seg. connecting to old patch on one end and new seg. on
C     other end integrate singular component (9) of surface current 
C     only.
   15 IF ( J1.LT.I1 .OR. J1.GT.I2 ) GOTO 18
      IPCH = ICON1(J1)
      IF ( IPCH.LT.PATOFF ) GOTO 16
      IPCH = IPCH - PATOFF
      FSIGN = -1.0D0
      GOTO 17
   16 IPCH = ICON2(J1)
      IF ( IPCH.LT.PATOFF ) GOTO 18
      IPCH = IPCH - PATOFF
      FSIGN = 1.0D0
   17 IF ( IPCH.GT.M1 ) GOTO 18
      JS = LDP - IPCH
      IPGND = 1
      T1XJ = T1X(JS)
      T1YJ = T1Y(JS)
      T1ZJ = T1Z(JS)
      T2XJ = T2X(JS)
      T2YJ = T2Y(JS)
      T2ZJ = T2Z(JS)
      XJ = X(JS)
      YJ = Y(JS)
      ZJ = Z(JS)
      S = BI(JS)
      XI = X(J1)
      YI = Y(J1)
      ZI = Z(J1)
      CABI = CAB(J1)
      SABI = SAB(J1)
      SALPI = SALP(J1)
      CALL PCINT(XI,YI,ZI,CABI,SABI,SALPI,EMEL)
      PY = PI*SI(J1)*FSIGN
      PX = DSIN(PY)
      PY = DCOS(PY)
      EXC = EMEL(9)*FSIGN
      IL = JCO(JSNO)
      K = J1 - I1 + 1
      CW(K,IL) = CW(K,IL) + EXC*(AX(JSNO)+BX(JSNO)*PX+CX(JSNO)*PY)
 
   18 RETURN
 
      END
