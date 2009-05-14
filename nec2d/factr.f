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
 
      SUBROUTINE FACTR( NN , A , IIP , NDIM , LD, D , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine to factor a matrix into a unit lower triangular matrix  **
C* and an upper triangular matrix using the Gauss-Doolittle algorithm **
C* presented on pages 411-416 of A. Ralston: A first Course In        **
C* Numerical Analysis. Comments below refer to comments in Ralstons   **
C* text. The matrix is this case is transposed.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NN                                                        **
C* OUTPUT - A                                                         **
C* OUTPUT - IIP                                                       **
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
C* called by   FACGF   FACTRS  NETWK                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NN , NDIM , LD , DEBUG
      INTEGER IIP(NDIM)
      COMPLEX*16 A(NDIM,NDIM) , D(2*LD)

C     Local variables.
      INTEGER I , IFLG , J , JP1 , K , R , RM1 , RP1 , PJ , PR
      REAL*8 DMAX , ELMAG
      COMPLEX*16 ARJ


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FACTR'
      ENDIF
 
      IFLG = 0
      DO 9 R = 1 , NN

C     Step 1.

         DO 1 K = 1 , NN
            D(K) = A(R,K)
    1    CONTINUE

C     Steps 2 and 3

         RM1 = R - 1
         IF ( RM1.LT.1 ) GOTO 4
         DO 3 J = 1 , RM1
            PJ = IIP(J)
            ARJ = D(PJ)
            A(R,J) = ARJ
            D(PJ) = D(J)
            JP1 = J + 1
            DO 2 I = JP1 , NN
               D(I) = D(I) - A(J,I)*ARJ
    2       CONTINUE
    3    CONTINUE
    4    CONTINUE

C     Step 4

         DMAX = DBLE(D(R)*DCONJG(D(R)))
         IIP(R) = R
         RP1 = R + 1
         IF ( RP1.GT.NN ) GOTO 6
         DO 5 I = RP1 , NN
            ELMAG = DBLE(D(I)*DCONJG(D(I)))
            IF ( ELMAG.LT.DMAX ) GOTO 5
            DMAX = ELMAG
            IIP(R) = I
    5    CONTINUE
    6    CONTINUE
         IF ( DMAX.LT.1.D-10 ) IFLG = 1
         PR = IIP(R)
         A(R,R) = D(PR)
         D(PR) = D(R)

C     Step 5

         IF ( RP1.GT.NN ) GOTO 8
         ARJ = 1.0D0/A(R,R)
         DO 7 I = RP1 , NN
            A(R,I) = D(I)*ARJ
    7    CONTINUE
    8    CONTINUE
         IF ( IFLG.EQ.0 ) GOTO 9
         WRITE (CHRSLT,10) R , DMAX
         IFLG = 0
    9 CONTINUE
 
      RETURN
 
   10 FORMAT (' ','PIVOT(',I3,')=',1P,E16.8)
 
      END
