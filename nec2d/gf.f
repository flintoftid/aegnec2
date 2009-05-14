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

      SUBROUTINE GF( ZK , CO , SSI )

C*--------------------------------------------------------------------**
C*                                                                    **
C* GF computes the integrand exp(jkr)/(kr) for numerical integration. **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ZK                                                        **
C* OUTPUT - CO                                                        **
C* OUTPUT - SSI                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  IJX     RKB2    ZPK                                    **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   INTX                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 ZK , CO , SSI 

C     Local variables.
      REAL*8 RK , RKS , ZDK

C     Common storage.
      INCLUDE 'tmi.inc'
 

      ZDK = ZK - ZPK
      RK = DSQRT(RKB2+ZDK*ZDK)
      SSI = DSIN(RK)/RK
      IF ( IJX ) 1 , 2 , 1
    1 CO = DCOS(RK)/RK
      RETURN
    2 IF ( RK.LT.0.2D0 ) GOTO 3
      CO = (DCOS(RK)-1.0D0)/RK
      RETURN
    3 RKS = RK*RK
      CO = ((-1.388888888888888889D-3*RKS+4.16666666666666666667D-2)
     &     *RKS-0.5D0)*RK
 
      RETURN
 
      END
