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
 
      SUBROUTINE HINTG( XI , YI , ZI , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* HINTG computes the near H field of a patch current in free space   **
C* or over ground.                                                    **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XI                                                        **
C* INPUT  - YI                                                        **
C* INPUT  - ZI                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS                                            **
C* uses value  EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS     IPERF   KSYMP   S       T1XJ    T1YJ   **
C*             T1ZJ    T2XJ    T2YJ    T2ZJ    XJ      YJ      ZJ     **
C*             ZRATI                                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CMSS    NHFLD                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER DEBUG
      REAL*8 XI , YI , ZI

C     Local variables.
      INTEGER IIP
      REAL*8 CR , CTH , FPI , PX , PY , R , RFL , RK , RSQ , RX , RY , 
     &       RZ , SR , T1ZR , T2ZR , TP , XYMAG
      COMPLEX*16 GAM , F1X , F1Y , F1Z , F2X , F2Y , F2Z , RRV , RRH

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
      INCLUDE 'gnd.inc'

C     Data initialisation.
      DATA FPI /12.56637061435917295385057D0/
      DATA TP /6.283185307179586476925286D0/
 
 
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING HINTG'
      ENDIF

      RX = XI - XJ
      RY = YI - YJ
      RFL = -1.0D0
      EXK = (0.0D0,0.0D0)
      EYK = (0.0D0,0.0D0)
      EZK = (0.0D0,0.0D0)
      EXS = (0.0D0,0.0D0)
      EYS = (0.0D0,0.0D0)
      EZS = (0.0D0,0.0D0)
      DO 5 IIP = 1 , KSYMP
         RFL = -RFL
         RZ = ZI - ZJ*RFL
         RSQ = RX*RX + RY*RY + RZ*RZ
         IF ( RSQ.LT.1.D-20 ) GOTO 5
         R = DSQRT(RSQ)
         RK = TP*R
         CR = DCOS(RK)
         SR = DSIN(RK)
         GAM = -(DCMPLX(CR,-SR)+RK*DCMPLX(SR,CR))/(FPI*RSQ*R)*S
         EXC = GAM*RX
         EYC = GAM*RY
         EZC = GAM*RZ
         T1ZR = T1ZJ*RFL
         T2ZR = T2ZJ*RFL
         F1X = EYC*T1ZR - EZC*T1YJ
         F1Y = EZC*T1XJ - EXC*T1ZR
         F1Z = EXC*T1YJ - EYC*T1XJ
         F2X = EYC*T2ZR - EZC*T2YJ
         F2Y = EZC*T2XJ - EXC*T2ZR
         F2Z = EXC*T2YJ - EYC*T2XJ
         IF ( IIP.EQ.1 ) GOTO 4
         IF ( IPERF.NE.1 ) GOTO 1
         F1X = -F1X
         F1Y = -F1Y
         F1Z = -F1Z
         F2X = -F2X
         F2Y = -F2Y
         F2Z = -F2Z
         GOTO 4
    1    XYMAG = DSQRT(RX*RX+RY*RY)
         IF ( XYMAG.GT.1.D-6 ) GOTO 2
         PX = 0.0D0
         PY = 0.0D0
         CTH = 1.0D0
         RRV = (1.0D0,0.0D0)
         GOTO 3
    2    PX = -RY/XYMAG
         PY = RX/XYMAG
         CTH = RZ/R
         RRV = CDSQRT(1.0D0-ZRATI*ZRATI*(1.0D0-CTH*CTH))
    3    RRH = ZRATI*CTH
         RRH = (RRH-RRV)/(RRH+RRV)
         RRV = ZRATI*RRV
         RRV = -(CTH-RRV)/(CTH+RRV)
         GAM = (F1X*PX+F1Y*PY)*(RRV-RRH)
         F1X = F1X*RRH + GAM*PX
         F1Y = F1Y*RRH + GAM*PY
         F1Z = F1Z*RRH
         GAM = (F2X*PX+F2Y*PY)*(RRV-RRH)
         F2X = F2X*RRH + GAM*PX
         F2Y = F2Y*RRH + GAM*PY
         F2Z = F2Z*RRH
    4    EXK = EXK + F1X
         EYK = EYK + F1Y
         EZK = EZK + F1Z
         EXS = EXS + F2X
         EYS = EYS + F2Y
         EZS = EZS + F2Z
    5 CONTINUE
 
      RETURN
 
      END
