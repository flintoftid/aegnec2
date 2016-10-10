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

      SUBROUTINE BLCKIN( AR , NUNIT , IX1 , IX2 , NBLKS , NEOF , IFAIL ,
     &                   DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* BLCKIN controls the reading and writing of matrix blocks on files  **
C* for the out-of-core matrix solution.                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - AR                                                        **
C* INPUT  - NUNIT                                                     **
C* INPUT  - IX1                                                       **
C* INPUT  - IX2                                                       **
C* INPUT  - NBLKS                                                     **
C* INPUT  - NEOF                                                      **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FACIO   FACTRS  GFIL    GFOUT   LTSOLV  LUNSCR         **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NUNIT , IX1 , IX2 , NBLKS , NEOF , IFAIL , DEBUG
      COMPLEX*16 AR(1)

C     Local variables.
      INTEGER I , I1 , I2 , J


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING BLCKIN'
      ENDIF
      
      I1 = (IX1+1)/2
      I2 = (IX2+1)/2
      DO 2 I = 1 , NBLKS
         READ (NUNIT,END=3) (AR(J),J=I1,I2)
    2 CONTINUE
      RETURN
    3 WRITE (CHRSLT,4) NUNIT , NBLKS , NEOF
      IF ( NEOF.NE.777 ) THEN
         IFAIL=3
         RETURN
      ENDIF
C     IDF commented out one line - attempt to modify const arg.
C     NEOF = 0
 
      RETURN
 
    4 FORMAT ('  EOF ON UNIT',I3,'  NBLKS= ',I3,'  NEOF= ',I5)
 
      END
