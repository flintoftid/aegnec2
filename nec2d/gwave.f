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

      SUBROUTINE GWAVE( ERV , EZV , ERH , EZH , EPH )

C*--------------------------------------------------------------------**
C*                                                                    **
C* GWAVE computes the electric field, including ground wave, of a     **
C* current element over a ground plane using formulas of K.A. Norton  **
C* (Proc. IRE, Sept., 1937, pp.1203,1236.)                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - ERV                                                       **
C* OUTPUT - EZV                                                       **
C* OUTPUT - ERH                                                       **
C* OUTPUT - EZH                                                       **
C* OUTPUT - EPH                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  R1      R2      U       U2      XX1     XX2     ZMH    **
C*             ZPH                                                    **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       FBAR                                                   **
C* called by   GFLD    SFLDS                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      COMPLEX*16 FBAR
      EXTERNAL FBAR

C     Dummy arguments.
      COMPLEX*16 ERV , EZV , ERH , EZH , EPH

C     Local variables.
      REAL*8 CPP , CPP2 , CPPP , CPPP2 , SPP , SPP2 , SPPP , 
     &       SPPP2
      REAL*8 TPJX(2) , ECONX(2)
      COMPLEX*16 TPJ , RK1 , RK2 , TT1 , TT2 , T3 , T4 , P1 , RV , 
     &           OMR , W , F , Q1 , RH , V , G , XR1 , XR2 , X1 , XXX2 , 
     &           X3 , X4 , X5 , X6 , X7 , ECON

C     Common storage.
      INCLUDE 'gwav.inc'

C     Equivalent common storage.
      EQUIVALENCE (TPJ,TPJX) , (ECON,ECONX)
 
C     Data initialisation.
      DATA TPJX /0.0D0 , 6.283185307179586476925286D0/
      DATA ECONX /0.0D0 , -188.3651567308853277340992D0/
 

      SPPP = ZMH/R1
      SPPP2 = SPPP*SPPP
      CPPP2 = 1.0D0 - SPPP2
      IF ( CPPP2.LT.1.D-20 ) CPPP2 = 1.D-20
      CPPP = DSQRT(CPPP2)
      SPP = ZPH/R2
      SPP2 = SPP*SPP
      CPP2 = 1.0D0 - SPP2
      IF ( CPP2.LT.1.D-20 ) CPP2 = 1.D-20
      CPP = DSQRT(CPP2)
      RK1 = -TPJ*R1
      RK2 = -TPJ*R2
      TT1 = 1.0D0 - U2*CPP2
      TT2 = CDSQRT(TT1)
      T3 = (1.0D0-1.0D0/RK1)/RK1
      T4 = (1.0D0-1.0D0/RK2)/RK2
      P1 = RK2*U2*TT1/(2.0D0*CPP2)
      RV = (SPP-U*TT2)/(SPP+U*TT2)
      OMR = 1.0D0 - RV
      W = 1.0D0/OMR
      W = (4.0D0,0.0D0)*P1*W*W
      F = FBAR(W)
      Q1 = RK2*TT1/(2.0D0*U2*CPP2)
      RH = (TT2-U*SPP)/(TT2+U*SPP)
      V = 1.0D0/(1.0D0+RH)
      V = (4.0D0,0.0D0)*Q1*V*V
      G = FBAR(V)
      XR1 = XX1/R1
      XR2 = XX2/R2
      X1 = CPPP2*XR1
      XXX2 = RV*CPP2*XR2
      X3 = OMR*CPP2*F*XR2
      X4 = U*TT2*SPP*2.0D0*XR2/RK2
      X5 = XR1*T3*(1.0D0-3.0D0*SPPP2)
      X6 = XR2*T4*(1.0D0-3.0D0*SPP2)
      EZV = (X1+XXX2+X3-X4-X5-X6)*ECON
      X1 = SPPP*CPPP*XR1
      XXX2 = RV*SPP*CPP*XR2
      X3 = CPP*OMR*U*TT2*F*XR2
      X4 = SPP*CPP*OMR*XR2/RK2
      X5 = 3.0D0*SPPP*CPPP*T3*XR1
      X6 = CPP*U*TT2*OMR*XR2/RK2*0.5D0
      X7 = 3.0D0*SPP*CPP*T4*XR2
      ERV = -(X1+XXX2-X3+X4-X5+X6-X7)*ECON
      EZH = -(X1-XXX2+X3-X4-X5-X6+X7)*ECON
      X1 = SPPP2*XR1
      XXX2 = RV*SPP2*XR2
      X4 = U2*TT1*OMR*F*XR2
      X5 = T3*(1.0D0-3.0D0*CPPP2)*XR1
      X6 = T4*(1.0D0-3.0D0*CPP2)*(1.0D0-U2*(1.0D0+RV)-U2*OMR*F)*XR2
      X7 = U2*CPP2*OMR*(1.0D0-1.0D0/RK2)*
     &     (F*(U2*TT1-SPP2-1.0D0/RK2)+1.0D0/RK2)*XR2
      ERH = (X1-XXX2-X4-X5+X6+X7)*ECON
      X1 = XR1
      XXX2 = RH*XR2
      X3 = (RH+1.0D0)*G*XR2
      X4 = T3*XR1
      X5 = T4*(1.0D0-U2*(1.0D0+RV)-U2*OMR*F)*XR2
      X6 = 0.5D0*U2*OMR*(F*(U2*TT1-SPP2-1.0D0/RK2)+1.0D0/RK2)*XR2/RK2
      EPH = -(X1-XXX2+X3-X4+X5+X6)*ECON
 
      RETURN
 
      END
