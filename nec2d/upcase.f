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

      SUBROUTINE UPCASE( INTEXT , OUTTXT , LENGTH , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* UPCASE finds the length of INTEXT and converts it to upper case.   **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - INTEXT                                                    **
C* OUTPUT - OUTTXT                                                    **
C* OUTPUT - LENGTH                                                    **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   NEC2D   PARSIT                                         **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Dummy arguments.
      CHARACTER*(*) INTEXT , OUTTXT
      INTEGER LENGTH , DEBUG

C     Local variables.
      INTEGER I , J
 
      
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING UPCASE'
      ENDIF
         
      LENGTH = LEN(INTEXT)
      DO 3000 I = 1 , LENGTH
         J = ICHAR(INTEXT(I:I))
         IF ( J.GE.96 ) J = J - 32
         OUTTXT(I:I) = CHAR(J)
 3000 CONTINUE
 
      RETURN
 
      END
