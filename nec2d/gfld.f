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

      SUBROUTINE GFLD( RHO , PHI , RZ , ETH , EPI , ERD , UX , KKSYMP ,
     &                 LD , SALP , AIR , AII , BIR , BII , CIR , 
     &                 CII , CUR , N , M , X , Y , Z , SI , CAB , 
     &                 SAB , XS , YS , ZS , Z2S )

C*--------------------------------------------------------------------**
C*                                                                    **
C* GFLD computes the radiated E field at intermediate distances from  **
C* a radiating structure including the ground wave.                   **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - RHO                                                       **
C* INPUT  - PHI                                                       **
C* INPUT  - RZ                                                        **
C* OUTPUT - ETH                                                       **
C* OUTPUT - EPI                                                       **
C* OUTPUT - ERD                                                       **
C* INPUT  - UX                                                        **
C* INPUT  - KKSYMP                                                    **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - AIR                                                       **
C* INPUT  - AII                                                       **
C* INPUT  - BIR                                                       **
C* INPUT  - BII                                                       **
C* INPUT  - CIR                                                       **
C* INPUT  - CII                                                       **
C* PASSED - CUR                                                       **
C* INPUT  - N                                                         **
C* PASSED - M                                                         **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* INPUT  - CAB                                                       **
C* INPUT  - SAB                                                       **
C* PASSED - XS                                                        **
C* PASSED - YS                                                        **
C* PASSED - ZS                                                        **
C* PASSED - Z2S                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    R1      R2      U       U2      XX1     XX2     ZMH    **
C*             ZPH                                                    **
C* uses value  U                                                      **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       FFLD    GWAVE                                          **
C* called by   RDPAT                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL FFLD , GWAVE

C     Dummy arguments.
      INTEGER KKSYMP , LD , N , M 
      REAL*8 RHO , PHI , RZ
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       CAB(LD) , SAB(LD) ,XS(LD) , YS(LD) , ZS(LD) , Z2S(LD)
      COMPLEX*16 ETH , EPI , ERD , UX
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      INTEGER I , K
      REAL*8 A , ARG , BB , BOO , BOT , C , CALP , CBET , CPH , DX , 
     &       DY , DZ , EL , OMEGA , PHX , PHY , PI , R , RFL , RHP , 
     &       RHS , RHX , RHY , RI , RIX , RIY , RIZ , RNX , RNY , 
     &       RNZ , RR , RX , RXYZ , RY , SBET , SILL , SPH , THET , 
     &       THX , THY , THZ , TOO , TOP , TP
      COMPLEX*16 CIX , CIY , CIZ , EXA , ERV , EZV , ERH , EPH ,EZH , 
     &           EX , EY 

C     Common storage.
      INCLUDE 'gwav.inc'

C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA TP /6.283185307179586476925286D0/
 

      R = DSQRT(RHO*RHO+RZ*RZ)
      IF ( KKSYMP.EQ.1 ) GOTO 1
      IF ( CDABS(UX).GT.0.5D0 ) GOTO 1
      IF ( R.GT.1.0D5 ) GOTO 1
      GOTO 4

C     Computation of space wave only.

    1 IF ( RZ.LT.1.0D-20 ) GOTO 2
      THET = DATAN(RHO/RZ)
      GOTO 3
    2 THET = PI*0.5D0
    3 CALL FFLD(THET,PHI,ETH,EPI,LD,SALP,AIR,AII,BIR,BII,CIR,CII,
     &          CUR,N,M,X,Y,Z,SI,CAB,SAB,XS,YS,ZS,Z2S)
      ARG = -TP*R
      EXA = DCMPLX(DCOS(ARG),DSIN(ARG))/R
      ETH = ETH*EXA
      EPI = EPI*EXA
      ERD = (0.0D0,0.0D0)
      RETURN

C     Computation of space and ground waves.

    4 U = UX
      U2 = U*U
      PHX = -DSIN(PHI)
      PHY = DCOS(PHI)
      RX = RHO*PHY
      RY = -RHO*PHX
      CIX = (0.0D0,0.0D0)
      CIY = (0.0D0,0.0D0)
      CIZ = (0.0D0,0.0D0)

C     Summation of field from individual segments.

      DO 17 I = 1 , N
         DX = CAB(I)
         DY = SAB(I)
         DZ = SALP(I)
         RIX = RX - X(I)
         RIY = RY - Y(I)
         RHS = RIX*RIX + RIY*RIY
         RHP = DSQRT(RHS)
         IF ( RHP.LT.1.D-6 ) GOTO 5
         RHX = RIX/RHP
         RHY = RIY/RHP
         GOTO 6
    5    RHX = 1.0D0
         RHY = 0.0D0
    6    CALP = 1.0D0 - DZ*DZ
         IF ( CALP.LT.1.D-6 ) GOTO 7
         CALP = DSQRT(CALP)
         CBET = DX/CALP
         SBET = DY/CALP
         CPH = RHX*CBET + RHY*SBET
         SPH = RHY*CBET - RHX*SBET
         GOTO 8
    7    CPH = RHX
         SPH = RHY
    8    EL = PI*SI(I)
         RFL = -1.0D0

C     Integration of (current)*(phase factor) over segment and image for
C     constant, sine, and cosine current distributions.

         DO 16 K = 1 , 2
            RFL = -RFL
            RIZ = RZ - Z(I)*RFL
            RXYZ = DSQRT(RIX*RIX+RIY*RIY+RIZ*RIZ)
            RNX = RIX/RXYZ
            RNY = RIY/RXYZ
            RNZ = RIZ/RXYZ
            OMEGA = -(RNX*DX+RNY*DY+RNZ*DZ*RFL)
            SILL = OMEGA*EL
            TOP = EL + SILL
            BOT = EL - SILL
            IF ( DABS(OMEGA).LT.1.D-7 ) GOTO 9
            A = 2.0D0*DSIN(SILL)/OMEGA
            GOTO 10
    9       A = (2.0D0-OMEGA*OMEGA*EL*EL/3.0D0)*EL
   10       IF ( DABS(TOP).LT.1.D-7 ) GOTO 11
            TOO = DSIN(TOP)/TOP
            GOTO 12
   11       TOO = 1.0D0 - TOP*TOP/6.0D0
   12       IF ( DABS(BOT).LT.1.D-7 ) GOTO 13
            BOO = DSIN(BOT)/BOT
            GOTO 14
   13       BOO = 1.0D0 - BOT*BOT/6.0D0
   14       BB = EL*(BOO-TOO)
            C = EL*(BOO+TOO)
            RR = A*AIR(I) + BB*BII(I) + C*CIR(I)
            RI = A*AII(I) - BB*BIR(I) + C*CII(I)
            ARG = TP*(X(I)*RNX+Y(I)*RNY+Z(I)*RNZ*RFL)
            EXA = DCMPLX(DCOS(ARG),DSIN(ARG))*DCMPLX(RR,RI)/TP
            IF ( K.EQ.2 ) GOTO 15
            XX1 = EXA
            R1 = RXYZ
            ZMH = RIZ
            GOTO 16
   15       XX2 = EXA
            R2 = RXYZ
            ZPH = RIZ
   16    CONTINUE

C     Call subroutine to compute the field of segment including ground
C     wave.

         CALL GWAVE(ERV,EZV,ERH,EZH,EPH)
         ERH = ERH*CPH*CALP + ERV*DZ
         EPH = EPH*SPH*CALP
         EZH = EZH*CPH*CALP + EZV*DZ
         EX = ERH*RHX - EPH*RHY
         EY = ERH*RHY + EPH*RHX
         CIX = CIX + EX
         CIY = CIY + EY
         CIZ = CIZ + EZH
   17 CONTINUE
      ARG = -TP*R
      EXA = DCMPLX(DCOS(ARG),DSIN(ARG))
      CIX = CIX*EXA
      CIY = CIY*EXA
      CIZ = CIZ*EXA
      RNX = RX/R
      RNY = RY/R
      RNZ = RZ/R
      THX = RNZ*PHY
      THY = -RNZ*PHX
      THZ = -RHO/R
      ETH = CIX*THX + CIY*THY + CIZ*THZ
      EPI = CIX*PHX + CIY*PHY
      ERD = CIX*RNX + CIY*RNY + CIZ*RNZ
 
      RETURN
 
      END
