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

      SUBROUTINE BLCKOT ( AR , NUNIT , IX1 , IX2 , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* BLCKOT controls the reading and writing of matrix blocks on files  **
C* for the out-of-core matrix solution.                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - AR                                                        **
C* INPUT  - NUNIT                                                     **
C* INPUT  - IX1                                                       **
C* INPUT  - IX2                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CMSET   FACIO   FACTRS  GFIL    GFOUT   LUNSCR         **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Dummy arguments.
      INTEGER NUNIT , IX1 , IX2 , DEBUG
      COMPLEX*16 AR(1)

C     Local variables.
      INTEGER I1 , I2 , J

      
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING BLCKOT'
      ENDIF

      I1 = (IX1+1)/2
      I2 = (IX2+1)/2
      WRITE (NUNIT) (AR(J),J=I1,I2)

      RETURN

      END
