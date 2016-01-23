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

      SUBROUTINE LUNSCR( A , NROW , NOP , IX , IIP , IU2 , IU3 , IU4 ,
     &                   IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Unscrambles the lower triangular matrix of the factored out of     **
C* core matrix and determines appropiate ordering of unknowns.        **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - A                                                         **
C* INPUT  - NROW                                                      **
C* INPUT  - NOP                                                       **
C* OUTPUT - IX                                                        **
C* INPUT  - IIP                                                       **
C* INPUT  - IU2                                                       **
C* INPUT  - IU3                                                       **
C* INPUT  - IU4                                                       **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  NBLSYM  NPSYM                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN  BLCKOT                                         **
C* called by   FACGF   FACTRS                                         **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL BLCKIN , BLCKOT

C     Dummy arguments.
      INTEGER NROW , NOP , IU2 , IU3 , IU4 , IFAIL , DEBUG
      INTEGER IX(NROW) , IIP(NROW)
      COMPLEX*16 A(NROW,1)

C     Local variables.
      INTEGER I , I1 , I2 , IPI , IPK , IXBLK1 , IXT , J , J2 , K , 
     &        K1 , KA , KK , NB1 , NM1
      COMPLEX*16 TEMP

C     Common storage.
      INCLUDE 'matpar.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING LUNSCR'
      ENDIF
 
      I1 = 1
      I2 = 2*NPSYM*NROW
      NM1 = NROW - 1
      REWIND IU2
      REWIND IU3
      REWIND IU4
      DO 9 KK = 1 , NOP
         KA = (KK-1)*NROW
         DO 4 IXBLK1 = 1 , NBLSYM
            CALL BLCKIN(A,IU2,I1,I2,1,121,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            K1 = (IXBLK1-1)*NPSYM + 2
            IF ( NM1.LT.K1 ) GOTO 3
            J2 = 0
            DO 2 K = K1 , NM1
               IF ( J2.LT.NPSYM ) J2 = J2 + 1
               IPK = IIP(K+KA)
               DO 1 J = 1 , J2
                  TEMP = A(K,J)
                  A(K,J) = A(IPK,J)
                  A(IPK,J) = TEMP
    1          CONTINUE
    2       CONTINUE
    3       CONTINUE
            CALL BLCKOT(A,IU3,I1,I2,DEBUG)
    4    CONTINUE
         DO 5 IXBLK1 = 1 , NBLSYM
            BACKSPACE IU3
            IF ( IXBLK1.NE.1 ) BACKSPACE IU3
            CALL BLCKIN(A,IU3,I1,I2,1,123,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            CALL BLCKOT(A,IU4,I1,I2,DEBUG)
    5    CONTINUE
         DO 6 I = 1 , NROW
            IX(I+KA) = I
    6    CONTINUE
         DO 7 I = 1 , NROW
            IPI = IIP(I+KA)
            IXT = IX(I+KA)
            IX(I+KA) = IX(IPI+KA)
            IX(IPI+KA) = IXT
    7    CONTINUE
         IF ( NOP.EQ.1 ) GOTO 9
         NB1 = NBLSYM - 1
C     Skip NB1 logical records forward.
         DO 8 IXBLK1 = 1 , NB1
            CALL BLCKIN(A,IU3,I1,I2,1,125,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
    8    CONTINUE
    9 CONTINUE
      REWIND IU2
      REWIND IU3
      REWIND IU4
 
      RETURN
 
      END
