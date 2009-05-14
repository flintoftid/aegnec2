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

      SUBROUTINE LFACTR( A , NROW , IX1 , IX2 , IIP , LD , D , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* LFACTR performs Gauss-Doolittle manipulations on the two blocks of **
C* the transposed matrix in core storage. The Gauss-Doolittle         **
C* algorithm is presented on pages 411-416 of A. Ralston: A First     **
C* Course In numerical Analysis. Comments below refer to comments in  **
C* Ralstons text.                                                     **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - A                                                         **
C* INPUT  - NROW                                                      **
C* INPUT  - IX1                                                       **
C* INPUT  - IX2                                                       **
C* OUTPUT - IIP                                                       **
C* INPUT  - LD                                                        **
C* OUTPUT - D                                                         **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  NBLSYM  NLSYM   NPSYM                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FACIO                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NROW , IX1 , IX2 , LD , DEBUG
      INTEGER IIP(NROW)
      COMPLEX*16 A(NROW,1) , D(2*LD)

C     Local variables.
      LOGICAL L1 , L2 , L3
      INTEGER I , IFLG , IXJ , J , J1 , J2 , J2P1 , J2P2 , JP1 , 
     &        K , R , RR1 , RR2 , PJ , PR
      REAL*8 DMAX , ELMAG
      COMPLEX*16 AJR

C     Common storage.
      INCLUDE 'matpar.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING LFACTR'
      ENDIF
 
      IFLG = 0

C     Initialize RR1,RR2,J1,J2.

      L1 = IX1.EQ.1 .AND. IX2.EQ.2
      L2 = (IX2-1).EQ.IX1
      L3 = IX2.EQ.NBLSYM
      IF ( L1 ) GOTO 1
      GOTO 2
    1 RR1 = 1
      RR2 = 2*NPSYM
      J1 = 1
      J2 = -1
      GOTO 5
    2 RR1 = NPSYM + 1
      RR2 = 2*NPSYM
      J1 = (IX1-1)*NPSYM + 1
      IF ( L2 ) GOTO 3
      GOTO 4
    3 J2 = J1 + NPSYM - 2
      GOTO 5
    4 J2 = J1 + NPSYM - 1
    5 IF ( L3 ) RR2 = NPSYM + NLSYM
      DO 16 R = RR1 , RR2

C     Step 1.

         DO 6 K = J1 , NROW
            D(K) = A(K,R)
    6    CONTINUE

C     Steps 2 and 3.

         IF ( L1 .OR. L2 ) J2 = J2 + 1
         IF ( J1.GT.J2 ) GOTO 9
         IXJ = 0
         DO 8 J = J1 , J2
            IXJ = IXJ + 1
            PJ = IIP(J)
            AJR = D(PJ)
            A(J,R) = AJR
            D(PJ) = D(J)
            JP1 = J + 1
            DO 7 I = JP1 , NROW
               D(I) = D(I) - A(I,IXJ)*AJR
    7       CONTINUE
    8    CONTINUE
    9    CONTINUE

C     Step 4.

         J2P1 = J2 + 1
         IF ( L1 .OR. L2 ) GOTO 11
         IF ( NROW.LT.J2P1 ) GOTO 16
         DO 10 I = J2P1 , NROW
            A(I,R) = D(I)
   10    CONTINUE
         GOTO 16
   11    DMAX = DBLE(D(J2P1)*DCONJG(D(J2P1)))
         IIP(J2P1) = J2P1
         J2P2 = J2 + 2
         IF ( J2P2.GT.NROW ) GOTO 13
         DO 12 I = J2P2 , NROW
            ELMAG = DBLE(D(I)*DCONJG(D(I)))
            IF ( ELMAG.LT.DMAX ) GOTO 12
            DMAX = ELMAG
            IIP(J2P1) = I
   12    CONTINUE
   13    CONTINUE
         IF ( DMAX.LT.1.D-10 ) IFLG = 1
         PR = IIP(J2P1)
         A(J2P1,R) = D(PR)
         D(PR) = D(J2P1)

C     Step 5.

         IF ( J2P2.GT.NROW ) GOTO 15
         AJR = 1.0D0/A(J2P1,R)
         DO 14 I = J2P2 , NROW
            A(I,R) = D(I)*AJR
   14    CONTINUE
   15    CONTINUE
         IF ( IFLG.EQ.0 ) GOTO 16
         WRITE (CHRSLT,17) J2 , DMAX
         IFLG = 0
   16 CONTINUE
 
      RETURN
 
   17 FORMAT (' ','PIVOT(',I3,')=',1P,E16.8)
 
      END
