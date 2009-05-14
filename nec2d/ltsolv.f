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

      SUBROUTINE LTSOLV( A , NROW , IX , BB , NNEQ , NRH , IFL1 , IFL2 ,
     &                   LD , D , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* LTSOLV solves the matrix eq. Y(R)*LU(T)=BB(R) where (R) denotes row**
C* vector and LU(T) denotes the LU decomposition of the transpose of  **
C* the original coefficient matrix. The LU(T) decomposition is        **
C* stored on tape 5 in blocks in ascending order and on file 3 in     **
C* blocks of descending order.                                        **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - A                                                         **
C* INPUT  - NROW                                                      **
C* INPUT  - IX                                                        **
C* OUTPUT - BB                                                        **
C* INPUT  - NNEQ                                                      **
C* INPUT  - NRH                                                       **
C* PASSED - IFL1                                                      **
C* PASSED - IFL2                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - D                                                         **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  NBLSYM  NLSYM   NPSYM                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN                                                 **
C* called by   SOLGF   SOLVES                                         **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL BLCKIN

C     Dummy arguments.
      INTEGER NROW , NNEQ , NRH , IFL1 , IFL2 , LD , IFAIL , DEBUG
      INTEGER IX(NNEQ)
      COMPLEX*16 A(NROW,NROW) , BB(NNEQ,NRH) , D(2*LD)

C     Local variables.
      INTEGER I , I2 , IC , IXBLK1 , IXI , J , JM1 , JP1 , JST , K , 
     &        K2 , KP 
      COMPLEX*16 SUM

C     Common storage.
      INCLUDE 'matpar.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING LTSOLV'
      ENDIF
 
C     Forward substitution.

      I2 = 2*NPSYM*NROW
      DO 5 IXBLK1 = 1 , NBLSYM
         CALL BLCKIN(A,IFL1,1,I2,1,121,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         K2 = NPSYM
         IF ( IXBLK1.EQ.NBLSYM ) K2 = NLSYM
         JST = (IXBLK1-1)*NPSYM
         DO 4 IC = 1 , NRH
            J = JST
            DO 3 K = 1 , K2
               JM1 = J
               J = J + 1
               SUM = DCMPLX(0.0D0,0.0D0)
               IF ( JM1.LT.1 ) GOTO 2
               DO 1 I = 1 , JM1
                  SUM = SUM + A(I,K)*BB(I,IC)
    1          CONTINUE
    2          BB(J,IC) = (BB(J,IC)-SUM)/A(J,K)
    3       CONTINUE
    4    CONTINUE
    5 CONTINUE

C     Backward substitution.

      JST = NROW + 1
      DO 9 IXBLK1 = 1 , NBLSYM
         CALL BLCKIN(A,IFL2,1,I2,1,122,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         K2 = NPSYM
         IF ( IXBLK1.EQ.1 ) K2 = NLSYM
         DO 8 IC = 1 , NRH
            KP = K2 + 1
            J = JST
            DO 7 K = 1 , K2
               KP = KP - 1
               JP1 = J
               J = J - 1
               SUM = DCMPLX(0.0D0,0.0D0)
               IF ( NROW.LT.JP1 ) GOTO 7
               DO 6 I = JP1 , NROW
                  SUM = SUM + A(I,KP)*BB(I,IC)
    6          CONTINUE
               BB(J,IC) = BB(J,IC) - SUM
    7       CONTINUE
    8    CONTINUE
         JST = JST - K2
    9 CONTINUE

C     Unscramble solution.

      DO 12 IC = 1 , NRH
         DO 10 I = 1 , NROW
            IXI = IX(I)
            D(IXI) = BB(I,IC)
   10    CONTINUE
         DO 11 I = 1 , NROW
            BB(I,IC) = D(I)
   11    CONTINUE
   12 CONTINUE
 
      RETURN
 
      END
