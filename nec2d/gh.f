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
 
      SUBROUTINE GH( ZK , HR , HI )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Integrand for H field of a wire.                                   **  
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ZK                                                        **
C* OUTPUT - HR                                                        **
C* OUTPUT - HI                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  RHKS    ZPK2                                           **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   HFK                                                    **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 ZK , HR , HI

C     Local variables.
      REAL*8 CKR , R , RR2 , RR3 , RS , SKR

C     Common storage.
      INCLUDE 'tmh.inc'
 

      RS = ZK - ZPK2
      RS = RHKS + RS*RS
      R = DSQRT(RS)
      CKR = DCOS(R)
      SKR = DSIN(R)
      RR2 = 1.0D0/RS
      RR3 = RR2/R
      HR = SKR*RR2 + CKR*RR3
      HI = CKR*RR2 - SKR*RR3
 
      RETURN
 
      END
