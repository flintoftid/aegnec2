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
 
      SUBROUTINE TEST( F1R , F2R , TR , F1I , F2I , TI , DMIN )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Test for convergence in numerical integration for Romberg variable **
C* interval width integration routines.                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - F1R                                                       **
C* INPUT  - F2R                                                       **
C* OUTPUT - TR                                                        **
C* INPUT  - F1I                                                       **
C* INPUT  - F2I                                                       **
C* OUTPUT - TI                                                        **
C* INPUT  - DMIN                                                      **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   HFK     INTX    ROM2                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 F1R , F2R , TR , F1I , F2I , TI , DMIN

C***  Local variables.
      REAL*8 DEN
 

      DEN = DABS(F2R)
      TR = DABS(F2I)
      IF ( DEN.LT.TR ) DEN = TR
      IF ( DEN.LT.DMIN ) DEN = DMIN
      IF ( DEN.LT.1.0D-37 ) GOTO 1
      TR = DABS((F1R-F2R)/DEN)
      TI = DABS((F1I-F2I)/DEN)
      RETURN
    1 TR = 0.0D0
      TI = 0.0D0
 
      RETURN
 
      END
