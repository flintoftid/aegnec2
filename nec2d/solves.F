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

#include <config.h>
 
C*--------------------------------------------------------------------**

      SUBROUTINE SOLVES( A , IIP , BB , NNEQ , NRH , NNP , NN , MMP , 
     &                   MM , IFL1 , IFL2 , LD , D , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine SOLVES, for symmetric structures, handles the           **
C* transformation of the right hand side vector and solution of the   **
C* matrix equation.                                                   **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - A                                                         **
C* PASSED - IIP                                                       **
C* OUTPUT - BB                                                        **
C* INPUT  - NNEQ                                                      **
C* INPUT  - NRH                                                       **
C* INPUT  - NNP                                                       **
C* INPUT  - NN                                                        **
C* INPUT  - MMP                                                       **
C* INPUT  - MM                                                        **
C* INPUT  - IFL1                                                      **
C* INPUT  - IFL2                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - D                                                         **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ICASE   SSX                                            **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       LTSOLV  SOLVE   ZGETRS                                 **
C* called by   FACGF   SOLGF                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL LTSOLV , SOLVE
#ifdef USE_BLAS
      EXTERNAL ZGETRS
#endif

C     Dummy arguments.
      INTEGER NNEQ , NRH , NNP , NN , MMP , MM , IFL1 , IFL2 , LD ,
     &        IFAIL , DEBUG
      INTEGER IIP(1)
      COMPLEX*16 A(1) , BB(NNEQ,NRH) , D(2*LD)

C     Local variables.
      INTEGER I , IA , IB , IC , J , K , KK , NOP , NNPEQ , NROW
      REAL*8 FNOP , FNORM
      COMPLEX*16 SUM
#ifdef USE_BLAS
      INTEGER INFO
#endif

C     Common storage.
      INCLUDE 'smat.inc'
      INCLUDE 'matpar.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SOLVES'
      ENDIF
 
      NNPEQ = NNP + 2*MMP
      NOP = NNEQ/NNPEQ
      FNOP = NOP
      FNORM = 1.0D0/FNOP
      NROW = NNEQ
      IF ( ICASE.GT.3 ) NROW = NNPEQ
      IF ( NOP.EQ.1 ) GOTO 13
      DO 12 IC = 1 , NRH
         IF ( NN.EQ.0 .OR. MM.EQ.0 ) GOTO 6
         DO 1 I = 1 , NNEQ
            D(I) = BB(I,IC)
    1    CONTINUE
         KK = 2*MMP
         IA = NNP
         IB = NN
         J = NNP
         DO 5 K = 1 , NOP
            IF ( K.EQ.1 ) GOTO 3
            DO 2 I = 1 , NNP
               IA = IA + 1
               J = J + 1
               BB(J,IC) = D(IA)
    2       CONTINUE
            IF ( K.EQ.NOP ) GOTO 5
    3       DO 4 I = 1 , KK
               IB = IB + 1
               J = J + 1
               BB(J,IC) = D(IB)
    4       CONTINUE
    5    CONTINUE

C     Transform matrix eq. RHS vector according to symmetry modes.

    6    DO 11 I = 1 , NNPEQ
            DO 7 K = 1 , NOP
               IA = I + (K-1)*NNPEQ
               D(K) = BB(IA,IC)
    7       CONTINUE
            SUM = D(1)
            DO 8 K = 2 , NOP
               SUM = SUM + D(K)
    8       CONTINUE
            BB(I,IC) = SUM*FNORM
            DO 10 K = 2 , NOP
               IA = I + (K-1)*NNPEQ
               SUM = D(1)
               DO 9 J = 2 , NOP
                  SUM = SUM + D(J)*DCONJG(SSX(K,J))
    9          CONTINUE
               BB(IA,IC) = SUM*FNORM
   10       CONTINUE
   11    CONTINUE
   12 CONTINUE
   13 IF ( ICASE.LT.3 ) GOTO 14
      REWIND IFL1
      REWIND IFL2

C     Solve each mode equation.

   14 DO 18 KK = 1 , NOP
         IA = (KK-1)*NNPEQ + 1
         IB = IA
         IF ( ICASE.NE.4 ) GOTO 15
         I = NNPEQ*NNPEQ
         READ (IFL1) (A(J),J=1,I)
         IB = 1
   15    IF ( ICASE.EQ.3 .OR. ICASE.EQ.5 ) GOTO 17

#ifdef USE_BLAS
         IF ( ICASE.EQ.4 ) THEN
            DO 16 IC = 1 , NRH
               CALL SOLVE(NNPEQ,A(IB),IIP(IA),BB(IA,IC),NROW,LD,D,
     &                    DEBUG)
 16         CONTINUE
         ELSE
            DO 916 IC = 1 , NRH
               CALL ZGETRS('T',NNPEQ,1,A(IB),NROW,IIP(IA),BB(IA,IC),
     &                     NROW,INFO)
 916        CONTINUE
         ENDIF
#else 
         DO 16 IC = 1 , NRH
            CALL SOLVE(NNPEQ,A(IB),IIP(IA),BB(IA,IC),NROW,LD,D,DEBUG)
   16    CONTINUE
#endif

         GOTO 18
   17    CALL LTSOLV(A,NNPEQ,IIP(IA),BB(IA,1),NNEQ,NRH,IFL1,IFL2,LD,
     &               D,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
   18 CONTINUE
      IF ( NOP.EQ.1 ) RETURN

C     Inverse transform the mode solutions.

      DO 29 IC = 1 , NRH
         DO 23 I = 1 , NNPEQ
            DO 19 K = 1 , NOP
               IA = I + (K-1)*NNPEQ
               D(K) = BB(IA,IC)
   19       CONTINUE
            SUM = D(1)
            DO 20 K = 2 , NOP
               SUM = SUM + D(K)
   20       CONTINUE
            BB(I,IC) = SUM
            DO 22 K = 2 , NOP
               IA = I + (K-1)*NNPEQ
               SUM = D(1)
               DO 21 J = 2 , NOP
                  SUM = SUM + D(J)*SSX(K,J)
   21          CONTINUE
               BB(IA,IC) = SUM
   22       CONTINUE
   23    CONTINUE
         IF ( NN.EQ.0 .OR. MM.EQ.0 ) GOTO 29
         DO 24 I = 1 , NNEQ
            D(I) = BB(I,IC)
   24    CONTINUE
         KK = 2*MMP
         IA = NNP
         IB = NN
         J = NNP
         DO 28 K = 1 , NOP
            IF ( K.EQ.1 ) GOTO 26
            DO 25 I = 1 , NNP
               IA = IA + 1
               J = J + 1
               BB(IA,IC) = D(J)
   25       CONTINUE
            IF ( K.EQ.NOP ) GOTO 28
   26       DO 27 I = 1 , KK
               IB = IB + 1
               J = J + 1
               BB(IB,IC) = D(J)
   27       CONTINUE
   28    CONTINUE
   29 CONTINUE
 
      RETURN
 
      END
