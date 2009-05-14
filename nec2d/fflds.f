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

      SUBROUTINE FFLDS( ROX , ROY , ROZ , SCUR , EX , EY ,EZ , LD ,
     &                  M , XS , YS , ZS , Z2S )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Calculates the xyz components of the electric field due to         **
C* surface currents.                                                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ROX                                                       **
C* INPUT  - ROY                                                       **
C* INPUT  - ROZ                                                       **
C* INPUT  - SCUR                                                      **
C* OUTPUT - EX                                                        **
C* OUTPUT - EY                                                        **
C* OUTPUT - EZ                                                        **
C* INPUT  - LD                                                        **
C* INPUT  - M                                                         **
C* INPUT  - XS                                                        **
C* INPUT  - YS                                                        **
C* INPUT  - ZS                                                        **
C* INPUT  - Z2S                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FFLD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER LD , M 
      REAL*8 ROX , ROY , ROZ 
      REAL*8 XS(LD) , YS(LD) , ZS(LD) , Z2S(LD)
      COMPLEX*16 SCUR(1) , EX , EY , EZ

C     Local variables.
      INTEGER I , J , K
      REAL*8 ARG , TPI 
      REAL*8 CONSX(2)
      COMPLEX*16 CT , CONS 

C     Equivalent local storage.
      EQUIVALENCE (CONS,CONSX)

C     Data initialisation
      DATA TPI /6.283185307179586476925286D0/
      DATA CONSX /0.0D0 , 188.3651567308853277340992D0/
 

      EX = (0.0D0,0.0D0)
      EY = (0.0D0,0.0D0)
      EZ = (0.0D0,0.0D0)
      I = LD + 1
      DO 1 J = 1 , M
         I = I - 1
         ARG = TPI*(ROX*XS(I)+ROY*YS(I)+ROZ*ZS(I))
         CT = DCMPLX(DCOS(ARG)*Z2S(I),DSIN(ARG)*Z2S(I))
         K = 3*J
         EX = EX + SCUR(K-2)*CT
         EY = EY + SCUR(K-1)*CT
         EZ = EZ + SCUR(K)*CT
    1 CONTINUE
      CT = ROX*EX + ROY*EY + ROZ*EZ
      EX = CONS*(CT*ROX-EX)
      EY = CONS*(CT*ROY-EY)
      EZ = CONS*(CT*ROZ-EZ)
 
      RETURN
 
      END
