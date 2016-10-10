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

      SUBROUTINE SFLDS( T , E )

C*--------------------------------------------------------------------**
C*                                                                    **
C* SFLDX returns the field due to ground for a current element on     **
C* the source segment at T relative to the segment center.            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - T                                                         **
C* OUTPUT - E                                                         **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    R1      R2      XX1     XX2     ZMH     ZPH            **
C* uses value  CABJ    FRATI   ISNOR   R2      S       SABJ    SALPJ  **
C*             SN      XJ      XO      XSN     XX2     YJ      YO     **
C*             YSN     ZJ      ZO      ZPH                            **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       GWAVE   INTRP                                          **
C* called by   EFLD    ROM2                                           **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL GWAVE , INTRP

C     Dummy arguments.
      REAL*8 T
      COMPLEX*16 E(9)

C     Local variables.
      REAL*8 CPH , PHX , PHY , PI , POT , R2S , RHO , RHS , RHX , 
     &       RHY , RK , SFAC , SPH , THET , TP , XT , YT , ZPHS , ZT
      COMPLEX*16 ERV , EZV , ERH , EZH , EPH , ER , ET , HRV , 
     &           HZV , HRH

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'incom.inc'
      INCLUDE 'gwav.inc'
      INCLUDE 'gnd.inc'

C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA TP /6.283185307179586476925286D0/
      DATA POT /1.570796326794896619231321D0/
 

      XT = XJ + T*CABJ
      YT = YJ + T*SABJ
      ZT = ZJ + T*SALPJ
      RHX = XO - XT
      RHY = YO - YT
      RHS = RHX*RHX + RHY*RHY
      RHO = DSQRT(RHS)
      IF ( RHO.GT.0. ) GOTO 1
      RHX = 1.0D0
      RHY = 0.0D0
      PHX = 0.0D0
      PHY = 1.0D0
      GOTO 2
    1 RHX = RHX/RHO
      RHY = RHY/RHO
      PHX = -RHY
      PHY = RHX
    2 CPH = RHX*XSN + RHY*YSN
      SPH = RHY*XSN - RHX*YSN
      IF ( DABS(CPH).LT.1.D-10 ) CPH = 0.0D0
      IF ( DABS(SPH).LT.1.D-10 ) SPH = 0.0D0
      ZPH = ZO + ZT
      ZPHS = ZPH*ZPH
      R2S = RHS + ZPHS
      R2 = DSQRT(R2S)
      RK = R2*TP
      XX2 = DCMPLX(DCOS(RK),-DSIN(RK))
      IF ( ISNOR.EQ.1 ) GOTO 3

C     Use norton approximation for field due to ground.  Current is
C     lumped at segment center with current moment for constant, sine,
C     or cosine distribution.

      ZMH = 1.0D0
      R1 = 1.0D0
      XX1 = 0.0D0
      CALL GWAVE(ERV,EZV,ERH,EZH,EPH)
      ET = -(0.0D0,4.77134D0)*FRATI*XX2/(R2S*R2)
      ER = 2.0D0*ET*DCMPLX(1.0D0,RK)
      ET = ET*DCMPLX(1.0D0-RK*RK,RK)
      HRV = (ER+ET)*RHO*ZPH/R2S
      HZV = (ZPHS*ER-RHS*ET)/R2S
      HRH = (RHS*ER-ZPHS*ET)/R2S
      ERV = ERV - HRV
      EZV = EZV - HZV
      ERH = ERH + HRH
      EZH = EZH + HRV
      EPH = EPH + ET
      ERV = ERV*SALPJ
      EZV = EZV*SALPJ
      ERH = ERH*SN*CPH
      EZH = EZH*SN*CPH
      EPH = EPH*SN*SPH
      ERH = ERV + ERH
      E(1) = (ERH*RHX+EPH*PHX)*S
      E(2) = (ERH*RHY+EPH*PHY)*S
      E(3) = (EZV+EZH)*S
      E(4) = 0.0D0
      E(5) = 0.0D0
      E(6) = 0.0D0
      SFAC = PI*S
      SFAC = DSIN(SFAC)/SFAC
      E(7) = E(1)*SFAC
      E(8) = E(2)*SFAC
      E(9) = E(3)*SFAC
      RETURN

C     Interpolate in Sommerfeld field tables.

    3 IF ( RHO.LT.1.0D-12 ) GOTO 4
      THET = DATAN(ZPH/RHO)
      GOTO 5
    4 THET = POT
    5 CALL INTRP(R2,THET,ERV,EZV,ERH,EPH)
C     Combine vertical and horizontal components and convert to X,Y,Z
C     Components.  Multiply by EXP(-JKR)/R.
      XX2 = XX2/R2
      SFAC = SN*CPH
      ERH = XX2*(SALPJ*ERV+SFAC*ERH)
      EZH = XX2*(SALPJ*EZV-SFAC*ERV)
      EPH = SN*SPH*XX2*EPH
C     X,Y,Z fields for constant current.
      E(1) = ERH*RHX + EPH*PHX
      E(2) = ERH*RHY + EPH*PHY
      E(3) = EZH
      RK = TP*T
C     X,Y,Z fields for sine current.
      SFAC = DSIN(RK)
      E(4) = E(1)*SFAC
      E(5) = E(2)*SFAC
      E(6) = E(3)*SFAC
C     X,Y,Z fields for cosine current.
      SFAC = DCOS(RK)
      E(7) = E(1)*SFAC
      E(8) = E(2)*SFAC
      E(9) = E(3)*SFAC
 
      RETURN
 
      END
