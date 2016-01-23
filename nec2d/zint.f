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

      COMPLEX*16 FUNCTION ZINT( SIGL , RMU , ROLAM , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* ZINT computes the internal impedance of a circular wire.           **
C* Fix from NEC81 to fix permeability.                                **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - SIGL                                                      **
C* INPUT  - RMU                                                       **
C* INPUT  - ROLAM                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   LOAD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER DEBUG
      REAL*8  SIGL , RMU , ROLAM

C     Local variables.
      REAL*8 BEI , BER , CMOTP , DD , PI , POT , SSS , TP , TPCMU , 
     &       XX , YY
      REAL*8 FJX(2) , CNX(2) , CCN(28)
      COMPLEX*16 TH , PH , F , G , FJ , CN , BR1 , BR2 , CC1 , CC2 , 
     &           CC3 , CC4 , CC5 , CC6 , CC7 , CC8 , CC9 , CC10 , 
     &           CC11 , CC12 , CC13 , CC14

C     Equivalent storage.
      EQUIVALENCE (FJ,FJX) , (CN,CNX) , (CC1,CCN(1)) , (CC2,CCN(3)) , 
     &            (CC3,CCN(5)) , (CC4,CCN(7)) , (CC5,CCN(9)) , 
     &            (CC6,CCN(11)) , (CC7,CCN(13)) , (CC8,CCN(15)) , 
     &            (CC9,CCN(17)) , (CC10,CCN(19)) , (CC11,CCN(21)) , 
     &            (CC12,CCN(23)) , (CC13,CCN(25)) , (CC14,CCN(27))
 
C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA POT /1.570796326794896619231321D0/
      DATA TP /6.283185307179586476925286D0/
      DATA TPCMU /2367.066370312157358387101D0/
      DATA CMOTP /59.9584916D0/
      DATA FJX /0.0D0 , 1.0D0/
      DATA CNX /0.7071067811865475244008443D0 , 
     &          0.7071067811865475244008443D0/
      DATA CCN /6.0D-7 , 1.9D-6 , -3.4D-6 , 5.1D-6 , -2.52D-5 , 0.0D0 , 
     &          -9.06D-5 , -9.01D-5 , 0.0D0 , -9.765D-4 , 0.0110486D0 , 
     &          -0.0110485D0 , 0.0D0 , -0.3926991D0 , 1.6D-6 , -3.2D-6 , 
     &          1.17D-5 , -2.4D-6 , 3.46D-5 , 3.38D-5 , 5.0D-7 , 
     &          2.452D-4 , -1.3813D-3 , 1.3811D-3 , -6.25001D-2 , 
     &          -1.0D-7 , 0.7071068D0 , 0.7071068D0/
  
      TH(DD) = (((((CC1*DD+CC2)*DD+CC3)*DD+CC4)*DD+CC5)*DD+CC6)*DD + CC7
      PH(DD) = (((((CC8*DD+CC9)*DD+CC10)*DD+CC11)*DD+CC12)*DD+CC13)*DD 
     &       + CC14
      F(DD) = DSQRT(POT/DD)*CDEXP(-CN*DD+TH(-8.0D0/XX))
      G(DD) = CDEXP(CN*DD+TH(8.0D0/XX))/DSQRT(TP*DD)


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING ZINT'
      ENDIF

C     IDF - fix from NEC81 for permeability.
C     XX=DSQRT(TPCMU*SIGL)*ROLAM
      XX = DSQRT(TPCMU*RMU*SIGL)*ROLAM
C     End of fix.
      IF ( XX.GT.110.0D0 ) GOTO 2
      IF ( XX.GT.8.0D0 ) GOTO 1
      YY = XX/8.0D0
      YY = YY*YY
      SSS = YY*YY
      BER = ((((((-9.01D-6*SSS+1.22552D-3)*SSS-.08349609D+0)*SSS
     &    + 2.6419140D+0)*SSS-32.363456D+0)*SSS+113.77778D+0)*SSS
     &    - 64.0D0)*SSS + 1.0D0
      BEI = ((((((1.1346D-4*SSS-0.01103667D0)*SSS+0.52185615D0)*SSS
     &    - 10.567658D0)*SSS+72.817777D0)*SSS-113.77778D0)
     &    * SSS+16.0D0)*YY
      BR1 = DCMPLX(BER,BEI)
      BER = (((((((-3.94D-6*SSS+4.5957D-4)*SSS-0.02609253D0)*SSS
     &    + 0.66047849D0)*SSS-6.0681481D0)*SSS+14.222222D0)*SSS
     &    - 4.0D0)*YY)*XX
      BEI = ((((((4.609D-5*SSS-3.79386D-3)*SSS+0.14677204D0)*SSS
     &    - 2.3116751D0)*SSS+11.377778D0)*SSS-10.666667D0)*SSS
     &    + 0.5D0)*XX
      BR2 = DCMPLX(BER,BEI)
      BR1 = BR1/BR2
      GOTO 3
    1 BR2 = FJ*F(XX)/PI
      BR1 = G(XX) + BR2
      BR2 = G(XX)*PH(8.0D0/XX) - BR2*PH(-8.0D0/XX)
      BR1 = BR1/BR2
      GOTO 3
    2 BR1 = DCMPLX(.70710678D+0,-.70710678D+0)
C     IDF - fix from NEC81 for permeability.
C3    ZINT=FJ*SQRT(CMOTP/SIGL)*BR1/ROLAM
    3 ZINT = FJ*DSQRT(CMOTP*RMU/SIGL)*BR1/ROLAM
C     End of fix.
 
      RETURN
 
      END
