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

      SUBROUTINE FBLOCK( NROW , NCOL , IMAX , IRNGF , IIPSYM , IFAIL ,
     &                   DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FBLOCK sets parameters for out-of-core solution for the primary    **
C* matrix (A).                                                        **
C*                                                                    **C
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NROW                                                      **
C* INPUT  - NCOL                                                      **
C* INPUT  - IMAX                                                      **
C* INPUT  - IRNGF                                                     **
C* INPUT  - IIPSYM                                                    **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ICASE   IMAT    NBLOKS  NBLSYM  NLAST   NLSYM   NPBLK  **
C*             NPSYM   SSX                                            **
C* uses value  NBLOKS  NBLSYM  NLAST   NLSYM   NPBLK   NPSYM   SSX    **
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
      INTEGER NROW , NCOL , IMAX , IRNGF , IIPSYM , IFAIL , DEBUG

C     Local variables.
      INTEGER I , IMX1 , J , K , KA , KK , NOP
      REAL*8 ARG , PHAZ , TP
      COMPLEX*16 DETER

C     Common storage.
      INCLUDE 'matpar.inc'
      INCLUDE 'smat.inc'
 
C     Data initialisation.
      DATA TP /6.283185307179586476925286D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FBLOCK'
      ENDIF
 
      IMX1 = IMAX - IRNGF
      IF ( NROW*NCOL.GT.IMX1 ) GOTO 2
      NBLOKS = 1
      NPBLK = NROW
      NLAST = NROW
      IMAT = NROW*NCOL
      IF ( NROW.NE.NCOL ) GOTO 1
      ICASE = 1
      RETURN
    1 ICASE = 2
      GOTO 5
    2 CONTINUE

C     Following 4 lines from NEC81.
C     OPEN (CHTMP1,FORM='UNFORMATTED')
C     OPEN (CHTMP2,FORM='UNFORMATTED')
C     OPEN (CHTMP3,FORM='UNFORMATTED')
C     OPEN (CHTMP4,FORM='UNFORMATTED')

      IF ( NROW.NE.NCOL ) GOTO 3
      ICASE = 3
      NPBLK = IMAX/(2*NCOL)
      NPSYM = IMX1/NCOL
      IF ( NPSYM.LT.NPBLK ) NPBLK = NPSYM
      IF ( NPBLK.LT.1 ) GOTO 14
      NBLOKS = (NROW-1)/NPBLK
      NLAST = NROW - NBLOKS*NPBLK
      NBLOKS = NBLOKS + 1
      NBLSYM = NBLOKS
      NPSYM = NPBLK
      NLSYM = NLAST
      IMAT = NPBLK*NCOL
      WRITE (CHRSLT,16) NBLOKS , NPBLK , NLAST
      GOTO 13
    3 NPBLK = IMAX/NCOL
      IF ( NPBLK.LT.1 ) GOTO 14
      IF ( NPBLK.GT.NROW ) NPBLK = NROW
      NBLOKS = (NROW-1)/NPBLK
      NLAST = NROW - NBLOKS*NPBLK
      NBLOKS = NBLOKS + 1
      WRITE (CHRSLT,16) NBLOKS , NPBLK , NLAST
      IF ( NROW*NROW.GT.IMX1 ) GOTO 4
      ICASE = 4
      NBLSYM = 1
      NPSYM = NROW
      NLSYM = NROW
      IMAT = NROW*NROW
      WRITE (CHRSLT,17)
      GOTO 5
    4 ICASE = 5
      NPSYM = IMAX/(2*NROW)
      NBLSYM = IMX1/NROW
      IF ( NBLSYM.LT.NPSYM ) NPSYM = NBLSYM
      IF ( NPSYM.LT.1 ) GOTO 14
      NBLSYM = (NROW-1)/NPSYM
      NLSYM = NROW - NBLSYM*NPSYM
      NBLSYM = NBLSYM + 1
      WRITE (CHRSLT,18) NBLSYM , NPSYM , NLSYM
      IMAT = NPSYM*NROW
    5 NOP = NCOL/NROW
      IF ( NOP*NROW.NE.NCOL ) GOTO 15
      IF ( IIPSYM.GT.0 ) GOTO 8

C     Set up SSX matrix for rotational symmetry.

      PHAZ = TP/NOP
      DO 7 I = 2 , NOP
         DO 6 J = I , NOP
            ARG = PHAZ*DBLE(I-1)*DBLE(J-1)
            SSX(I,J) = DCMPLX(DCOS(ARG),DSIN(ARG))
            SSX(J,I) = SSX(I,J)
    6    CONTINUE
    7 CONTINUE
      GOTO 13

C     Set up SSX matrix for plane symmetry.

    8 KK = 1
      SSX(1,1) = DCMPLX(1.0D0,0.0D0)
      IF ( (NOP.EQ.2) .OR. (NOP.EQ.4) .OR. (NOP.EQ.8) ) GOTO 9
      IFAIL=15
      RETURN
    9 KA = NOP/2
      IF ( NOP.EQ.8 ) KA = 3
      DO 12 K = 1 , KA
         DO 11 I = 1 , KK
            DO 10 J = 1 , KK
               DETER = SSX(I,J)
               SSX(I,J+KK) = DETER
               SSX(I+KK,J+KK) = -DETER
               SSX(I+KK,J) = DETER
   10       CONTINUE
   11    CONTINUE
         KK = KK*2
   12 CONTINUE
   13 RETURN
   14 WRITE (CHRSLT,19) NROW , NCOL
      IFAIL=16
      RETURN
   15 WRITE (CHRSLT,20) NROW , NCOL
      IFAIL=17
      RETURN
 
   16 FORMAT (//' MATRIX FILE STORAGE -  NO. BLOCKS=',I5,
     &        ' COLUMNS PER BLOCK=',I5,' COLUMNS IN LAST BLOCK=',I5)
   17 FORMAT (' SUBMATRICIES FIT IN CORE')
   18 FORMAT (' SUBMATRIX PARTITIONING -  NO. BLOCKS=',I5,
     &        ' COLUMNS PER BLOCK=',I5,' COLUMNS IN LAST BLOCK=',I5)
   19 FORMAT (' ERROR - INSUFFICIENT STORAGE FOR MATRIX',2I5)
   20 FORMAT (' SYMMETRY ERROR - NROW,NCOL=',2I5)
 
      END
