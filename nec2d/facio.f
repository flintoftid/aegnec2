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

      SUBROUTINE FACIO( A , NROW , NOP , IIP , IU1 , IU2 , IU3 , IU4 , 
     &                  LD , D , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FACIO reads and writes matrix blocks needed for out-of-core LU     **
C* decomposition.                                                     **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - A                                                         **
C* INPUT  - NROW                                                      **
C* INPUT  - NOP                                                       **
C* PASSED - IIP                                                       **
C* INPUT  - IU1                                                       **
C* INPUT  - IU2                                                       **
C* INPUT  - IU3                                                       **
C* INPUT  - IU4                                                       **
C* INPUT  - LD                                                        **
C* PASSED - D                                                         **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  NBLSYM  NPSYM                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN  BLCKOT  LFACTR  TIMER                          **
C* called by   FACGF   FACTRS                                         **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL BLCKIN , BLCKOT , LFACTR , TIMER

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NROW , NOP , IU1 , IU2 , IU3 , IU4 , LD , IFAIL , DEBUG
      INTEGER IIP(NROW)
      COMPLEX*16 A(NROW,1) , D(2*LD)

C     Local variables.
      INTEGER I1 , I2 , I3 , I4 , IFILE3 , IFILE4 , IT , IXBLK1 , 
     &        IXBLK2 , IXBP , KA , KK , NBM
      REAL*8 TT1 , TT2 , TIME

C     Common storage.
      INCLUDE 'matpar.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FACIO'
      ENDIF
  
      IT = 2*NPSYM*NROW
      NBM = NBLSYM - 1
      I1 = 1
      I2 = IT
      I3 = I2 + 1
      I4 = 2*IT
      TIME = 0.0D0
      REWIND IU1
      REWIND IU2
      DO 3 KK = 1 , NOP
         KA = (KK-1)*NROW + 1
         IFILE3 = IU1
         IFILE4 = IU3
         DO 2 IXBLK1 = 1 , NBM
            REWIND IU3
            REWIND IU4
            CALL BLCKIN(A,IFILE3,I1,I2,1,17,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IXBP = IXBLK1 + 1
            DO 1 IXBLK2 = IXBP , NBLSYM
               CALL BLCKIN(A,IFILE3,I3,I4,1,18,IFAIL,DEBUG)
               IF(IFAIL.NE.0) RETURN
               CALL TIMER(TT1)
               CALL LFACTR(A,NROW,IXBLK1,IXBLK2,IIP(KA),LD,D,DEBUG)
               CALL TIMER(TT2)
               TIME = TIME + TT2 - TT1
               IF ( IXBLK2.EQ.IXBP ) CALL BLCKOT(A,IU2,I1,I2,DEBUG)
               IF ( IXBLK1.EQ.NBM .AND. IXBLK2.EQ.NBLSYM ) IFILE4 = IU2
               CALL BLCKOT(A,IFILE4,I3,I4,DEBUG)
    1       CONTINUE
            IFILE3 = IU3
            IFILE4 = IU4
            IF ( (IXBLK1/2)*2.NE.IXBLK1 ) GOTO 2
            IFILE3 = IU4
            IFILE4 = IU3
    2    CONTINUE
    3 CONTINUE
      REWIND IU1
      REWIND IU2
      REWIND IU3
      REWIND IU4
      WRITE (CHRSLT,4) TIME
 
      RETURN
 
    4 FORMAT (' CP TIME TAKEN FOR FACTORIZATION = ',1P,E12.5)
 
      END
