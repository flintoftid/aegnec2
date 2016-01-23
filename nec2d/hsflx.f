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
 
      SUBROUTINE HSFLX( SSS , RH , ZPX , HPK , HPS , HPC )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Calculates H field of sine cosine, and constant current of segment.**
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - SSS                                                       **
C* INPUT  - RH                                                        **
C* INPUT  - ZPX                                                       **
C* OUTPUT - HPK                                                       **
C* OUTPUT - HPS                                                       **
C* OUTPUT - HPC                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       HFK                                                    **
C* called by   HSFLD                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL HFK

C     Dummy arguments.
      REAL*8 SSS , RH , ZPX
      COMPLEX*16 HPK , HPS , HPC

C     Local variables.
      REAL*8 CDK , DH , DK , HKI , HKR , HSS , PI8 , RR1 , RR2 , 
     &       RH2 , RHZ , SDK , TP , Z1 , ZZ2 , ZP
      REAL*8 FJX(2) , FJKX(2)
      COMPLEX*16 FJ , FJK , EKR1 , EKR2 , TT1 , TT2 , CONS

C     Equivalent storage.
      EQUIVALENCE (FJ,FJX) , (FJK,FJKX)
 
C     Data initialisation.
      DATA TP /6.283185307179586476925286D0/
      DATA FJX /0.0D0 , 1.0D0/
      DATA FJKX /0.0D0 , -6.283185307179586476925286D0/
      DATA PI8 /25.13274122871834590770114D0/
 

      IF ( RH.LT.1.D-10 ) GOTO 6
      IF ( ZPX.LT.0. ) GOTO 1
      ZP = ZPX
      HSS = 1.0D0
      GOTO 2
    1 ZP = -ZPX
      HSS = -1.0D0
    2 DH = 0.5D0*SSS
      Z1 = ZP + DH
      ZZ2 = ZP - DH
      IF ( ZZ2.LT.1.D-7 ) GOTO 3
      RHZ = RH/ZZ2
      GOTO 4
    3 RHZ = 1.0D0
    4 DK = TP*DH
      CDK = DCOS(DK)
      SDK = DSIN(DK)
      CALL HFK(-DK,DK,RH*TP,ZP*TP,HKR,HKI)
      HPK = DCMPLX(HKR,HKI)
      IF ( RHZ.LT.1.D-3 ) GOTO 5
      RH2 = RH*RH
      RR1 = DSQRT(RH2+Z1*Z1)
      RR2 = DSQRT(RH2+ZZ2*ZZ2)
      EKR1 = CDEXP(FJK*RR1)
      EKR2 = CDEXP(FJK*RR2)
      TT1 = Z1*EKR1/RR1
      TT2 = ZZ2*EKR2/RR2
      HPS = (CDK*(EKR2-EKR1)-FJ*SDK*(TT2+TT1))*HSS
      HPC = -SDK*(EKR2+EKR1) - FJ*CDK*(TT2-TT1)
      CONS = -FJ/(2.0D0*TP*RH)
      HPS = CONS*HPS
      HPC = CONS*HPC
      RETURN
    5 EKR1 = DCMPLX(CDK,SDK)/(ZZ2*ZZ2)
      EKR2 = DCMPLX(CDK,-SDK)/(Z1*Z1)
      TT1 = TP*(1.0D0/Z1-1.0D0/ZZ2)
      TT2 = CDEXP(FJK*ZP)*RH/PI8
      HPS = TT2*(TT1+(EKR1+EKR2)*SDK)*HSS
      HPC = TT2*(-FJ*TT1+(EKR1-EKR2)*CDK)
      RETURN
    6 HPS = (0.0D0,0.0D0)
      HPC = (0.0D0,0.0D0)
      HPK = (0.0D0,0.0D0)
 
      RETURN
 
      END
