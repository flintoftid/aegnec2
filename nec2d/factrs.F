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

#include <config.h>
 
C*--------------------------------------------------------------------**

      SUBROUTINE FACTRS( NNP , NROW , A , IIP , IX , IU1 , IU2 , IU3 ,
     &                   IU4 , LD , D , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FACTRS, for symmetric structure, transforms submatricies to form   **
C* matricies of the symmetric modes and calls routine to factor       **
C* matricies. If no symmetry, the routine is called to factor the     **
C* complete matrix.                                                   **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NNP                                                       **
C* INPUT  - NROW                                                      **
C* OUTPUT - A                                                         **
C* PASSED - IIP                                                       **
C* PASSED - IX                                                        **
C* INPUT  - IU1                                                       **
C* INPUT  - IU2                                                       **
C* INPUT  - IU3                                                       **
C* PASSED - IU4                                                       **
C* INPUT  - LD                                                        **
C* PASSED - D                                                         **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ICASE   NBLOKS  NBLSYM  NLAST   NLSYM   NPBLK   NPSYM  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN  BLCKOT  FACIO   FACTR   LUNSCR  ZGETRF         **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL BLCKIN , BLCKOT , FACIO , FACTR , LUNSCR
#ifdef USE_BLAS
      EXTERNAL ZGETRF
#endif

C     Dummy arguments.
      INTEGER NNP , NROW , IU1 , IU2 , IU3 , IU4 , LD , IFAIL , DEBUG
      INTEGER IIP(NROW) , IX(NROW)
      COMPLEX*16 A(1) , D(2*LD)

C     Local variables.
      INTEGER I , I2 , ICOLDX , ICOLS , IR1 , IR2 , IRR1 , IRR2 ,  J , 
     &        J2 , K , KA , KK , L , NOP 
#ifdef USE_BLAS
      INTEGER INFO
#endif

C     Common storage.
      INCLUDE 'matpar.inc'
 

      IF ( DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FACTRS'
      ENDIF
 
      NOP = NROW/NNP
      IF ( ICASE.GT.2 ) GOTO 2
      DO 1 KK = 1 , NOP
         KA = (KK-1)*NNP + 1
#ifdef USE_BLAS
         CALL ZGETRF(NNP,NNP,A(KA),NROW,IIP(KA),INFO)
#else
         CALL FACTR(NNP,A(KA),IIP(KA),NROW,LD,D,DEBUG)
#endif
    1 CONTINUE
      RETURN
    2 IF ( ICASE.GT.3 ) GOTO 3

C     Factor submatricies, or factor complete matrix if no symmetry
C     exists.

      CALL FACIO(A,NROW,NOP,IX,IU1,IU2,IU3,IU4,LD,D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL LUNSCR(A,NROW,NOP,IIP,IX,IU2,IU3,IU4,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      RETURN

C     Rewrite the matrices by columns on tape 13.

    3 I2 = 2*NPBLK*NROW
      REWIND IU2
      DO 7 K = 1 , NOP
         REWIND IU1
         ICOLS = NPBLK
         IR2 = K*NNP
         IR1 = IR2 - NNP + 1
         DO 6 L = 1 , NBLOKS
            IF ( NBLOKS.EQ.1 .AND. K.GT.1 ) GOTO 4
            CALL BLCKIN(A,IU1,1,I2,1,602,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( L.EQ.NBLOKS ) ICOLS = NLAST
    4       IRR1 = IR1
            IRR2 = IR2
            DO 5 ICOLDX = 1 , ICOLS
               WRITE (IU2) (A(I),I=IRR1,IRR2)
               IRR1 = IRR1 + NROW
               IRR2 = IRR2 + NROW
    5       CONTINUE
    6    CONTINUE
    7 CONTINUE
      REWIND IU1
      REWIND IU2
      IF ( ICASE.EQ.5 ) GOTO 10
      REWIND IU3
      IRR1 = NNP*NNP
      DO 9 KK = 1 , NOP
         IR1 = 1 - NNP
         IR2 = 0
         DO 8 I = 1 , NNP
            IR1 = IR1 + NNP
            IR2 = IR2 + NNP
            READ (IU2) (A(J),J=IR1,IR2)
    8    CONTINUE
         KA = (KK-1)*NNP + 1
         CALL FACTR(NNP,A,IIP(KA),NNP,LD,D,DEBUG)
         WRITE (IU3) (A(I),I=1,IRR1)
    9 CONTINUE
      REWIND IU2
      REWIND IU3
      RETURN
   10 I2 = 2*NPSYM*NNP
      DO 13 KK = 1 , NOP
         J2 = NPSYM
         DO 12 L = 1 , NBLSYM
            IF ( L.EQ.NBLSYM ) J2 = NLSYM
            IR1 = 1 - NNP
            IR2 = 0
            DO 11 J = 1 , J2
               IR1 = IR1 + NNP
               IR2 = IR2 + NNP
               READ (IU2) (A(I),I=IR1,IR2)
   11       CONTINUE
            CALL BLCKOT(A,IU1,1,I2,DEBUG)
   12    CONTINUE
   13 CONTINUE
      REWIND IU1
      CALL FACIO(A,NNP,NOP,IX,IU1,IU2,IU3,IU4,LD,D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL LUNSCR(A,NNP,NOP,IIP,IX,IU2,IU3,IU4,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
 
      RETURN
 
      END
