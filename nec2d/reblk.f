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

      SUBROUTINE REBLK( BB , BBX , NB , NBX , N2C , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Reblock array B in NGF solution from blocks of rows on to blocks   **
C* columns.                                                           **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - BB                                                        **
C* OUTPUT - BBX                                                       **
C* INPUT  - NB                                                        **
C* INPUT  - NBX                                                       **
C* INPUT  - N2C                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  NBBL    NBBX    NLBL    NLBX    NPBL    NPBX           **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FACGF                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER NB , NBX , N2C , DEBUG
      COMPLEX*16 BB(NB,1) , BBX(NBX,1)

C     Parameter definiotns.
      INCLUDE 'nec2d.inc'

C     Local variables.
      INTEGER I , IB , IBX , IX , J , NIB , NIX , NPB , NPX

C     Common storage.
      INCLUDE 'matpar.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING REBLK'
      ENDIF

      REWIND CHTMP6
      NIB = 0
      NPB = NPBL
      DO 4 IB = 1 , NBBL
         IF ( IB.EQ.NBBL ) NPB = NLBL
         REWIND CHTMP4
         NIX = 0
         NPX = NPBX
         DO 3 IBX = 1 , NBBX
            IF ( IBX.EQ.NBBX ) NPX = NLBX
            READ (CHTMP4) ((BBX(I,J),I=1,NPX),J=1,N2C)
            DO 2 I = 1 , NPX
               IX = I + NIX
               DO 1 J = 1 , NPB
                  BB(IX,J) = BBX(I,J+NIB)
    1          CONTINUE
    2       CONTINUE
            NIX = NIX + NPBX
    3    CONTINUE
         WRITE (CHTMP6) ((BB(I,J),I=1,NB),J=1,NPB)
         NIB = NIB + NPBL
    4 CONTINUE
      REWIND CHTMP4
      REWIND CHTMP6
 
      RETURN
 
      END
