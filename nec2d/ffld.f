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
C*--------------------------------------------------------------------**

      SUBROUTINE FFLD( THET , PHI , ETH , EPH , LD , SALP , AIR , 
     &                 AII , BIR , BII , CIR , CII , CUR , N , M , X , 
     &                 Y , Z , SI , CAB , SAB , XS , YS , ZS , Z2S )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FFLD calculates the far zone radiated electric fields,             **
C* due to currents on wires and patches in free space or over ground. **
C* The factor exp(j*k*r)/(r/lamda) not included.                      **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - THET                                                      **
C* INPUT  - PHI                                                       **
C* OUTPUT - ETH                                                       **
C* OUTPUT - EPH                                                       **
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
C* INPUT  - M                                                         **
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
C* modifies    ** NOTHING **                                          **
C* uses value  CH      CL      IFAR    IPERF   KSYMP   SCRWL   T1     **
C*             T2      ZRATI   ZRATI2                                 **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       FFLDS                                                  **
C* called by   GFLD    RDPAT                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL FFLDS

C     Dummy arguments.
      INTEGER LD , N , M
      REAL*8 THET , PHI 
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       CAB(LD), SAB(LD) ,XS(LD) , YS(LD) , ZS(LD) , Z2S(LD)
      COMPLEX*16 ETH , EPH
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      INTEGER I , IIP , K
      REAL*8 A , ARG , BB , BOO , BOT , C , DD , DARG , DR , EL , ETA , 
     &       OMEGA , PHX , PHY , PI , RFL , RI , ROX , ROY , ROZ , 
     &       ROZS , RR , RRZ , SILL , THX , THY , THZ , TOO , 
     &       TOP , TP , TTHET
      REAL*8 CONSX(2)
      COMPLEX*16 CIX , CIY , CIZ , EXA , CONST , CCX , CCY , CCZ , 
     &           CDP , ZRSIN , RRV , RRH , RRV1 , RRH1 , RRV2 , RRH2 , 
     &           TIX , TIY , TIZ , ZSCRN , EX , EY , EZ , EGX , EGY , 
     &           EGZ

C     Common storage.
      INCLUDE 'gnd.inc'

C     Equivalent local storage.
      EQUIVALENCE (CONST,CONSX)

C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA TP /6.283185307179586476925286D0/
      DATA ETA /376.7303134617706554681984D0/
      DATA CONSX /0.0D0 , -29.9792458D0/
 

C     Following five lines added by IDF to suppress compiler warnings.
      DARG = 0.0D0
      TTHET = 0.0D0
      CIX = (0.0D0,0.0D0)
      CIY = (0.0D0,0.0D0)
      CIZ = (0.0D0,0.0D0)
      
C
      PHX = -DSIN(PHI)
      PHY = DCOS(PHI)
      ROZ = DCOS(THET)
      ROZS = ROZ
      THX = ROZ*PHY
      THY = -ROZ*PHX
      THZ = -DSIN(THET)
      ROX = -THZ*PHY
      ROY = THZ*PHX
      IF ( N.EQ.0 ) GOTO 20

C     Loop for structure image if any.

      DO 19 K = 1 , KSYMP

C     Calculation of reflection coeffecients.

         IF ( K.EQ.1 ) GOTO 4
         IF ( IPERF.NE.1 ) GOTO 1

C     For perfect ground.

         RRV = -(1.0D0,0.0D0)
         RRH = -(1.0D0,0.0D0)
         GOTO 2

C     For infinite planar ground.

    1    ZRSIN = CDSQRT(1.0D0-ZRATI*ZRATI*THZ*THZ)
         RRV = -(ROZ-ZRATI*ZRSIN)/(ROZ+ZRATI*ZRSIN)
         RRH = (ZRATI*ROZ-ZRSIN)/(ZRATI*ROZ+ZRSIN)
    2    IF ( IFAR.LE.1 ) GOTO 3

C     For the cliff problem, two reflction coefficients calculated.

         RRV1 = RRV
         RRH1 = RRH
         TTHET = DTAN(THET)
         IF ( IFAR.EQ.4 ) GOTO 3
         ZRSIN = CDSQRT(1.0D0-ZRATI2*ZRATI2*THZ*THZ)
         RRV2 = -(ROZ-ZRATI2*ZRSIN)/(ROZ+ZRATI2*ZRSIN)
         RRH2 = (ZRATI2*ROZ-ZRSIN)/(ZRATI2*ROZ+ZRSIN)
         DARG = -TP*2.0D0*CH*ROZ
    3    ROZ = -ROZ
         CCX = CIX
         CCY = CIY
         CCZ = CIZ
    4    CIX = (0.0D0,0.0D0)
         CIY = (0.0D0,0.0D0)
         CIZ = (0.0D0,0.0D0)

C     Loop over structure segments.

         DO 17 I = 1 , N
            OMEGA = -(ROX*CAB(I)+ROY*SAB(I)+ROZ*SALP(I))
            EL = PI*SI(I)
            SILL = OMEGA*EL
            TOP = EL + SILL
            BOT = EL - SILL
            IF ( DABS(OMEGA).LT.1.D-7 ) GOTO 5
            A = 2.0D0*DSIN(SILL)/OMEGA
            GOTO 6
    5       A = (2.0D0-OMEGA*OMEGA*EL*EL/3.0D0)*EL
    6       IF ( DABS(TOP).LT.1.D-7 ) GOTO 7
            TOO = DSIN(TOP)/TOP
            GOTO 8
    7       TOO = 1.0D0 - TOP*TOP/6.0D0
    8       IF ( DABS(BOT).LT.1.0D-7 ) GOTO 9
            BOO = DSIN(BOT)/BOT
            GOTO 10
    9       BOO = 1.0D0 - BOT*BOT/6.0D0
   10       BB = EL*(BOO-TOO)
            C = EL*(BOO+TOO)
            RR = A*AIR(I) + BB*BII(I) + C*CIR(I)
            RI = A*AII(I) - BB*BIR(I) + C*CII(I)
            ARG = TP*(X(I)*ROX+Y(I)*ROY+Z(I)*ROZ)
            IF ( K.EQ.2 .AND. IFAR.GE.2 ) GOTO 11
            EXA = DCMPLX(DCOS(ARG),DSIN(ARG))*DCMPLX(RR,RI)

C     Summation for far field integral.

            CIX = CIX + EXA*CAB(I)
            CIY = CIY + EXA*SAB(I)
            CIZ = CIZ + EXA*SALP(I)
            GOTO 17

C     Calculation of image contribution in cliff and ground screen
C     problems.

   11       DR = Z(I)*TTHET

C     Specular point distance.

            DD = DR*PHY + X(I)
            IF ( IFAR.EQ.2 ) GOTO 13
            DD = DSQRT(DD*DD+(Y(I)-DR*PHX)**2)
            IF ( IFAR.EQ.3 ) GOTO 13
            IF ( (SCRWL-DD).LT.0. ) GOTO 12

C     Radial wire ground screen reflection coefficient.

            DD = DD + T2
            ZSCRN = T1*DD*DLOG(DD/T2)
            ZSCRN = (ZSCRN*ZRATI)/(ETA*ZRATI+ZSCRN)
            ZRSIN = CDSQRT(1.0D0-ZSCRN*ZSCRN*THZ*THZ)
            RRV = (ROZ+ZSCRN*ZRSIN)/(-ROZ+ZSCRN*ZRSIN)
            RRH = (ZSCRN*ROZ+ZRSIN)/(ZSCRN*ROZ-ZRSIN)
            GOTO 16
   12       IF ( IFAR.EQ.4 ) GOTO 14
            IF ( IFAR.EQ.5 ) DD = DR*PHY + X(I)
   13       IF ( (CL-DD).LE.0. ) GOTO 15
   14       RRV = RRV1
            RRH = RRH1
            GOTO 16
   15       RRV = RRV2
            RRH = RRH2
            ARG = ARG + DARG
   16       EXA = DCMPLX(DCOS(ARG),DSIN(ARG))*DCMPLX(RR,RI)

C     Contribution of each image segment modified by reflection coef. ,
C     for cliff and ground screen problems.

            TIX = EXA*CAB(I)
            TIY = EXA*SAB(I)
            TIZ = EXA*SALP(I)
            CDP = (TIX*PHX+TIY*PHY)*(RRH-RRV)
            CIX = CIX + TIX*RRV + CDP*PHX
            CIY = CIY + TIY*RRV + CDP*PHY
            CIZ = CIZ - TIZ*RRV
   17    CONTINUE
         IF ( K.EQ.1 ) GOTO 19
         IF ( IFAR.GE.2 ) GOTO 18

C     Calculation of contribution of structure image for infinite 
C     ground.

         CDP = (CIX*PHX+CIY*PHY)*(RRH-RRV)
         CIX = CCX + CIX*RRV + CDP*PHX
         CIY = CCY + CIY*RRV + CDP*PHY
         CIZ = CCZ - CIZ*RRV
         GOTO 19
   18    CIX = CIX + CCX
         CIY = CIY + CCY
         CIZ = CIZ + CCZ
   19 CONTINUE
      IF ( M.GT.0 ) GOTO 21
      ETH = (CIX*THX+CIY*THY+CIZ*THZ)*CONST
      EPH = (CIX*PHX+CIY*PHY)*CONST
      RETURN
   20 CIX = (0.0D0,0.0D0)
      CIY = (0.0D0,0.0D0)
      CIZ = (0.0D0,0.0D0)
   21 ROZ = ROZS

C     Electric field components.

      RFL = -1.0D0
      DO 25 IIP = 1 , KSYMP
         RFL = -RFL
         RRZ = ROZ*RFL
         CALL FFLDS(ROX,ROY,RRZ,CUR(N+1),EGX,EGY,EGZ,LD,M,XS,YS,
     &              ZS,Z2S)
         IF ( IIP.EQ.2 ) GOTO 22
         EX = EGX
         EY = EGY
         EZ = EGZ
         GOTO 25
   22    IF ( IPERF.NE.1 ) GOTO 23
         EGX = -EGX
         EGY = -EGY
         EGZ = -EGZ
         GOTO 24
   23    RRV = CDSQRT(1.0D0-ZRATI*ZRATI*THZ*THZ)
         RRH = ZRATI*ROZ
         RRH = (RRH-RRV)/(RRH+RRV)
         RRV = ZRATI*RRV
         RRV = -(ROZ-RRV)/(ROZ+RRV)
         ETH = (EGX*PHX+EGY*PHY)*(RRH-RRV)
         EGX = EGX*RRV + ETH*PHX
         EGY = EGY*RRV + ETH*PHY
         EGZ = EGZ*RRV
   24    EX = EX + EGX
         EY = EY + EGY
         EZ = EZ - EGZ
   25 CONTINUE
      EX = EX + CIX*CONST
      EY = EY + CIY*CONST
      EZ = EZ + CIZ*CONST
      ETH = EX*THX + EY*THY + EZ*THZ
      EPH = EX*PHX + EY*PHY
 
      RETURN
 
      END
