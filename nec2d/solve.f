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

      SUBROUTINE SOLVE( NN , A , IIP , BB , NDIM , LD , D , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine to solve the matrix equation LU*X=BB where L is a unit  **
C* lower triangular matrix and U is an upper triangular matrix both   **
C* of which are stored in A. The RHS vector BB is input and the       **
C* solution is returned through vector BB. (Matrix transposed.)       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NN                                                        **
C* INPUT  - A                                                         **
C* INPUT  - IIP                                                       **
C* OUTPUT - BB                                                        **
C* INPUT  - NDIM                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - D                                                         **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   NETWK   SOLGF   SOLVES                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER NN , NDIM , LD , DEBUG
      INTEGER IIP(NDIM)
      COMPLEX*16 A(NDIM,NDIM) , BB(NDIM) , D(2*LD)

C     Local variables.
      INTEGER I , IP1 , J , K , PI
      COMPLEX*16 SUM


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SOLVE'
      ENDIF

C     Forward substitution.

      DO 3 I = 1 , NN
         PI = IIP(I)
         D(I) = BB(PI)
         BB(PI) = BB(I)
         IP1 = I + 1
         IF ( IP1.GT.NN ) GOTO 2
         DO 1 J = IP1 , NN
            BB(J) = BB(J) - A(I,J)*D(I)
    1    CONTINUE
    2    CONTINUE
    3 CONTINUE

C     Backward substitution.

      DO 6 K = 1 , NN
         I = NN - K + 1
         SUM = DCMPLX(0.0D0,0.0D0)
         IP1 = I + 1
         IF ( IP1.GT.NN ) GOTO 5
         DO 4 J = IP1 , NN
            SUM = SUM + A(J,I)*BB(J)
    4    CONTINUE
    5    CONTINUE
         BB(I) = (D(I)-SUM)/A(I,I)
    6 CONTINUE
 
      RETURN
 
      END
