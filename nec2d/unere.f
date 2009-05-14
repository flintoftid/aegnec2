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

      SUBROUTINE UNERE( XOB , YOB ,ZOB )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Calculates the electric field due to unit current in the t1 and t2 **
C* directions on a patch.                                             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XOB                                                       **
C* INPUT  - YOB                                                       **
C* INPUT  - ZOB                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXK     EXS     EYK     EYS     EZK     EZS            **
C* uses value  EXK     EXS     EYK     EYS     EZK     EZS     IPERF  **
C*             IPGND   S       T1XJ    T1YJ    T1ZJ    T2XJ    T2YJ   **
C*             T2ZJ    XJ      YJ      ZJ      ZRATI                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CMSW    NEFLD   PCINT                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 XOB , YOB , ZOB

C     Local variables.
      REAL*8 CONST , CTH , PX , PY , R , RR2 , RT , RX , RY , RZ , 
     &       T1ZR , T2ZR , TPI , TT1 , TT2 , XYMAG , ZR
      COMPLEX*16 ER , Q1 , Q2 , RRV , RRH , EDP

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
      INCLUDE 'gnd.inc'

C     Data initialisation.
       DATA TPI /6.283185307179586476925286D0/
       DATA CONST /4.771345159236942258888898D0/
 

C     CONST=ETA/(8.*PI**2)
      ZR = ZJ
      T1ZR = T1ZJ
      T2ZR = T2ZJ
      IF ( IPGND.NE.2 ) GOTO 1
      ZR = -ZR
      T1ZR = -T1ZR
      T2ZR = -T2ZR
    1 RX = XOB - XJ
      RY = YOB - YJ
      RZ = ZOB - ZR
      RR2 = RX*RX + RY*RY + RZ*RZ
      IF ( RR2.GT.1.D-20 ) GOTO 2
      EXK = (0.0D0,0.0D0)
      EYK = (0.0D0,0.0D0)
      EZK = (0.0D0,0.0D0)
      EXS = (0.0D0,0.0D0)
      EYS = (0.0D0,0.0D0)
      EZS = (0.0D0,0.0D0)
      RETURN
    2 R = DSQRT(RR2)
      TT1 = -TPI*R
      TT2 = TT1*TT1
      RT = RR2*R
      ER = DCMPLX(DSIN(TT1),-DCOS(TT1))*(CONST*S)
      Q1 = DCMPLX(TT2-1.0D0,TT1)*ER/RT
      Q2 = DCMPLX(3.0D0-TT2,-3.0D0*TT1)*ER/(RT*RR2)
      ER = Q2*(T1XJ*RX+T1YJ*RY+T1ZR*RZ)
      EXK = Q1*T1XJ + ER*RX
      EYK = Q1*T1YJ + ER*RY
      EZK = Q1*T1ZR + ER*RZ
      ER = Q2*(T2XJ*RX+T2YJ*RY+T2ZR*RZ)
      EXS = Q1*T2XJ + ER*RX
      EYS = Q1*T2YJ + ER*RY
      EZS = Q1*T2ZR + ER*RZ
      IF ( IPGND.EQ.1 ) GOTO 6
      IF ( IPERF.NE.1 ) GOTO 3
      EXK = -EXK
      EYK = -EYK
      EZK = -EZK
      EXS = -EXS
      EYS = -EYS
      EZS = -EZS
      GOTO 6
    3 XYMAG = DSQRT(RX*RX+RY*RY)
      IF ( XYMAG.GT.1.D-6 ) GOTO 4
      PX = 0.0D0
      PY = 0.0D0
      CTH = 1.0D0
      RRV = (1.0D0,0.0D0)
      GOTO 5
    4 PX = -RY/XYMAG
      PY = RX/XYMAG
      CTH = RZ/DSQRT(XYMAG*XYMAG+RZ*RZ)
      RRV = CDSQRT(1.0D0-ZRATI*ZRATI*(1.0D0-CTH*CTH))
    5 RRH = ZRATI*CTH
      RRH = (RRH-RRV)/(RRH+RRV)
      RRV = ZRATI*RRV
      RRV = -(CTH-RRV)/(CTH+RRV)
      EDP = (EXK*PX+EYK*PY)*(RRH-RRV)
      EXK = EXK*RRV + EDP*PX
      EYK = EYK*RRV + EDP*PY
      EZK = EZK*RRV
      EDP = (EXS*PX+EYS*PY)*(RRH-RRV)
      EXS = EXS*RRV + EDP*PX
      EYS = EYS*RRV + EDP*PY
      EZS = EZS*RRV
 
    6 RETURN
 
      END
