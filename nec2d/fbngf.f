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

      SUBROUTINE FBNGF( NNEQ , NNEQ2 , IB11 , IC11 , ID11 ,IX11 , 
     &                  IRESRV , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FBNGF sets the blocking parameters for the B, C, and D arrays for  **
C* out-of-core storage.                                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NNEQ                                                      **
C* INPUT  - NNEQ2                                                     **
C* OUTPUT - IB11                                                      **
C* OUTPUT - IC11                                                      **
C* OUTPUT - ID11                                                      **
C* OUTPUT - IX11                                                      **
C* INPUT  - IRESRV                                                    **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ICASX   NBBL    NBBX    NLBL    NLBX    NPBL    NPBX   **
C* uses value  ICASE   ICASX   IMAT    NBBL    NBBX    NLBL    NLBX   **
C*             NPBL    NPBX                                           **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NNEQ , NNEQ2 , IB11 , IC11 , ID11 , IX11 , IRESRV , 
     &        IFAIL , DEBUG
 
C     Local variables.
      INTEGER IR , IRESX , NBCD , NBLN , NDLN

C     Common storage.
      INCLUDE 'matpar.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FBNGF'
      ENDIF
 
      IRESX = IRESRV - IMAT
      NBLN = NNEQ*NNEQ2
      NDLN = NNEQ2*NNEQ2
      NBCD = 2*NBLN + NDLN
      IF ( NBCD.GT.IRESX ) GOTO 1
      ICASX = 1
      IB11 = IMAT + 1
      GOTO 2
    1 CONTINUE
C     Following 5 lines from NEC81.
C     OPEN (CHTMP1,FORM='UNFORMATTED')
C     OPEN (CHTMP2,FORM='UNFORMATTED')
C     OPEN (CHTMP3,FORM='UNFORMATTED')
C     OPEN (CHTMP4,FORM='UNFORMATTED')
C     OPEN (CHTMP5,FORM='UNFORMATTED')
C
      IF ( ICASE.LT.3 ) GOTO 3
      IF ( NBCD.GT.IRESRV .OR. NBLN.GT.IRESX ) GOTO 3
      ICASX = 2
      IB11 = 1
    2 NBBX = 1
      NPBX = NNEQ
      NLBX = NNEQ
      NBBL = 1
      NPBL = NNEQ2
      NLBL = NNEQ2
      GOTO 5
    3 IR = IRESRV
      IF ( ICASE.LT.3 ) IR = IRESX
      ICASX = 3
C     Following line from NEC81.
C     OPEN (CHTMP6,FORM='UNFORMATTED')
C
      IF ( NDLN.GT.IR ) ICASX = 4
      NBCD = 2*NNEQ + NNEQ2
      NPBL = IR/NBCD
      NLBL = IR/(2*NNEQ2)
      IF ( NLBL.LT.NPBL ) NPBL = NLBL
      IF ( ICASE.LT.3 ) GOTO 4
      NLBL = IRESX/NNEQ
      IF ( NLBL.LT.NPBL ) NPBL = NLBL
    4 IF ( NPBL.LT.1 ) GOTO 6
      NBBL = (NNEQ2-1)/NPBL
      NLBL = NNEQ2 - NBBL*NPBL
      NBBL = NBBL + 1
      NBLN = NNEQ*NPBL
      IR = IR - NBLN
      NPBX = IR/NNEQ2
      IF ( NPBX.GT.NNEQ ) NPBX = NNEQ
      NBBX = (NNEQ-1)/NPBX
      NLBX = NNEQ - NBBX*NPBX
      NBBX = NBBX + 1
      IB11 = 1
      IF ( ICASE.LT.3 ) IB11 = IMAT + 1
    5 IC11 = IB11 + NBLN
      ID11 = IC11 + NBLN
      IX11 = IMAT + 1
      WRITE (CHRSLT,11) NNEQ2
      IF ( ICASX.EQ.1 ) RETURN
      WRITE (CHRSLT,8) ICASX
      WRITE (CHRSLT,9) NBBX , NPBX , NLBX
      WRITE (CHRSLT,10) NBBL , NPBL , NLBL
      RETURN
    6 WRITE (CHRSLT,7) IRESRV , IMAT , NNEQ , NNEQ2
 
      IFAIL=18
      RETURN
 
    7 FORMAT (' ERROR - INSUFFICIENT STORAGE FOR INTERACTION MATRICIES',
     &        '  IRESRV,IMAT,NEQ,NEQ2 =',4I5)
    8 FORMAT (' FILE STORAGE FOR NEW MATRIX SECTIONS -  ICASX =',I2)
    9 FORMAT (' B FILLED BY ROWS -',15X,'NO. BLOCKS =',I3,3X,
     &        'ROWS PER BLOCK =',I3,3X,'ROWS IN LAST BLOCK =',I3)
   10 FORMAT (' B BY COLUMNS, C AND D BY ROWS -',2X,'NO. BLOCKS =',I3,
     &        4X,'R/C PER BLOCK =',I3,4X,'R/C IN LAST BLOCK =',I3)
   11 FORMAT (//,' N.G.F. - NUMBER OF NEW UNKNOWNS IS',I4)
 
      END
