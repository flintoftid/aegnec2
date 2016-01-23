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

      SUBROUTINE FACGF( A , BB , C , DD , BBX , IIP , IX , NNP , NN1 , 
     &                  MMP , MM1 , N1C , N2C , LD , D , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FACGF computes and factors DD-C(inv(A)BB) for the NGF solution.    **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - A                                                         **
C* OUTPUT - BB                                                        **
C* OUTPUT - C                                                         **
C* OUTPUT - DD                                                        **
C* OUTPUT - BBX                                                       **
C* PASSED - IIP                                                       **
C* PASSED - IX                                                        **
C* PASSED - NNP                                                       **
C* PASSED - NN1                                                       **
C* PASSED - MMP                                                       **
C* PASSED - MM1                                                       **
C* INPUT  - N1C                                                       **
C* INPUT  - N2C                                                       **
C* INPUT  - LD                                                        **
C* PASSED - D                                                         **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ICASE   NBLSYM  NLSYM   NPSYM                          **
C* passes arg  NPBX                                                   **
C* uses value  ICASE   ICASX   NBBL    NBLSYM  NLBL    NLSYM   NPBL   **
C*             NPBX    NPSYM                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       FACIO   FACTR   LUNSCR  REBLK   SOLVES                 **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL FACIO , FACTR , LUNSCR , REBLK , SOLVES

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NNP , NN1 , MMP , MM1 , N1C , N2C , LD , IFAIL , DEBUG
      INTEGER IIP(1) , IX(1)
      COMPLEX*16 A(1) , BB(N1C,1) , C(N1C,1) , DD(N2C,1) , BBX(N1C,1) ,
     &           D(2*LD)

C     Local variables.
      INTEGER I , IB , IBFL , IC , ICASS , II , J , K , N1CP , NBLSYS , 
     &        NIC , NLSYS , NPB , NPC , NPSYS
      COMPLEX*16 SUM

C     Common storage.
      INCLUDE 'matpar.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FACGF'
      ENDIF
 
      IF ( N2C.EQ.0 ) RETURN
      IBFL = CHTMP4
      IF ( ICASX.LT.3 ) GOTO 1
C     Convert BB from blocks of rows on CHTMP4 to blocks of columns on 
C     CHTMP6.
      CALL REBLK(BB,C,N1C,NPBX,N2C,DEBUG)
      IBFL = CHTMP6
    1 NPB = NPBL
      IF ( ICASX.EQ.2 ) REWIND CHTMP4
C     Compute INV(A)BB and write on CHTMP4.
      DO 2 IB = 1 , NBBL
         IF ( IB.EQ.NBBL ) NPB = NLBL
         IF ( ICASX.GT.1 ) READ (IBFL) ((BBX(I,J),I=1,N1C),J=1,NPB)
         CALL SOLVES(A,IIP,BBX,N1C,NPB,NNP,NN1,MMP,MM1,CHTMP3,CHTMP3,
     &               LD,D,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF ( ICASX.EQ.2 ) REWIND CHTMP4
         IF ( ICASX.GT.1 ) WRITE (CHTMP4) ((BBX(I,J),I=1,N1C),J=1,NPB)
    2 CONTINUE
      IF ( ICASX.EQ.1 ) GOTO 3
      REWIND CHTMP1
      REWIND CHTMP2
      REWIND CHTMP5
      REWIND IBFL
    3 NPC = NPBL
C     Compute DD-C(INV(A)BB) and write on CHTMP11.
      DO 9 IC = 1 , NBBL
         IF ( IC.EQ.NBBL ) NPC = NLBL
         IF ( ICASX.EQ.1 ) GOTO 4
         READ (CHTMP5) ((C(I,J),I=1,N1C),J=1,NPC)
         READ (CHTMP2) ((DD(I,J),I=1,N2C),J=1,NPC)
         REWIND CHTMP4
    4    NPB = NPBL
         NIC = 0
         DO 8 IB = 1 , NBBL
            IF ( IB.EQ.NBBL ) NPB = NLBL
            IF ( ICASX.GT.1 ) READ (CHTMP4) ((BB(I,J),I=1,N1C),J=1,NPB)
            DO 7 I = 1 , NPB
               II = I + NIC
               DO 6 J = 1 , NPC
                  SUM = DCMPLX(0.0D0,0.0D0)
                  DO 5 K = 1 , N1C
                     SUM = SUM + BB(K,I)*C(K,J)
    5             CONTINUE
                  DD(II,J) = DD(II,J) - SUM
    6          CONTINUE
    7       CONTINUE
            NIC = NIC + NPBL
    8    CONTINUE
         IF ( ICASX.GT.1 ) WRITE (CHTMP1) ((DD(I,J),I=1,N2C),J=1,NPBL)
    9 CONTINUE
      IF ( ICASX.EQ.1 ) GOTO 10
      REWIND CHTMP1
      REWIND CHTMP2
      REWIND CHTMP4
      REWIND CHTMP5
   10 N1CP = N1C + 1
C     Factor DD-C(INV(A)BB).
      IF ( ICASX.GT.1 ) GOTO 11
      CALL FACTR(N2C,DD,IIP(N1CP),N2C,LD,D,DEBUG)
      GOTO 14
   11 IF ( ICASX.EQ.4 ) GOTO 13
      NPB = NPBL
      IC = 0
      DO 12 IB = 1 , NBBL
         IF ( IB.EQ.NBBL ) NPB = NLBL
         II = IC + 1
         IC = IC + N2C*NPB
         READ (CHTMP1) (BB(I,1),I=II,IC)
   12 CONTINUE
      REWIND CHTMP1
      CALL FACTR(N2C,BB,IIP(N1CP),N2C,LD,D,DEBUG)
      NIC = N2C*N2C
      WRITE (CHTMP1) (BB(I,1),I=1,NIC)
      REWIND CHTMP1
      GOTO 14
   13 NBLSYS = NBLSYM
      NPSYS = NPSYM
      NLSYS = NLSYM
      ICASS = ICASE
      NBLSYM = NBBL
      NPSYM = NPBL
      NLSYM = NLBL
      ICASE = 3
      CALL FACIO(BB,N2C,1,IX(N1CP),CHTMP1,CHTMP2,CHTMP6,CHTMP1,LD,
     &           D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL LUNSCR(BB,N2C,1,IIP(N1CP),IX(N1CP),CHTMP2,CHTMP1,CHTMP6,
     &            IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN 
      NBLSYM = NBLSYS
      NPSYM = NPSYS
      NLSYM = NLSYS
      ICASE = ICASS
 
   14 RETURN
 
      END
