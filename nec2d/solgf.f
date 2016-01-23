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

      SUBROUTINE SOLGF( ZA , ZB , ZC , ZD , XY , IIP , NNP , NN1 , NN , 
     &                  MMP , MM1 , MM , N1C , N2C , N2CZ , LD , D ,
     &                  NSMAXX , NSCON , NPCON , ISCON , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Solve for current basis functions in NGF procedure.                **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - ZA                                                        **
C* OUTPUT - ZB                                                        **
C* OUTPUT - ZC                                                        **
C* PASSED - ZD                                                        **
C* OUTPUT - XY                                                        **
C* PASSED - IIP                                                       **
C* PASSED - NNP                                                       **
C* INPUT  - NN1                                                       **
C* INPUT  - NN                                                        **
C* PASSED - MMP                                                       **
C* INPUT  - MM1                                                       **
C* PASSED - MM                                                        **
C* INPUT  - N1C                                                       **
C* INPUT  - N2C                                                       **
C* INPUT  - N2CZ                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - D                                                         **
C* INPUT  - NSMAXX                                                    **
C* INPUT  - NSCON                                                     **
C* INPUT  - NPCON                                                     **
C* INPUT  - ISCON                                                     **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ICASE   NBLSYM  NLSYM   NPSYM                          **
C* uses value  ICASE   ICASX   NBBL    NBLSYM  NLBL    NLSYM   NPBL   **
C*             NPSYM                                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       LTSOLV  SOLVE   SOLVES                                 **
C* called by   NETWK                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL LTSOLV , SOLVE , SOLVES

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NNP , NN , NN1 , MMP , MM , MM1 , N1C , N2C , N2CZ , 
     &        LD , NSMAXX , NSCON , NPCON , IFAIL , DEBUG
      INTEGER IIP(1) , ISCON(NSMAXX)
      COMPLEX*16 ZA(1) , ZB(N1C,1) , ZC(N1C,1) , ZD(N2CZ,1) , XY(1) ,
     &           D(2*LD)

C     Local variables.
      INTEGER I , ICASS , IFL , II , J , JJ , JP , NN2 , NBLSYS , 
     &        NNEQ , NEQS , NI , NLSYS , NPB , NPM , NPSYS
      COMPLEX*16 SUM

C     Common storage.
      INCLUDE 'matpar.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SOLGF'
      ENDIF
 
      IFL = 14
      IF ( ICASX.GT.0 ) IFL = 13
      IF ( N2C.GT.0 ) GOTO 1
C     Normal solution.  Not N.G.F.
      CALL SOLVES(ZA,IIP,XY,N1C,1,NNP,NN,MMP,MM,13,IFL,LD,D,IFAIL,
     &            DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 22
    1 IF ( NN1.EQ.NN .OR. MM1.EQ.0 ) GOTO 5
C     Reorder excitation array.
      NN2 = NN1 + 1
      JJ = NN + 1
      NPM = NN + 2*MM1
      DO 2 I = NN2 , NPM
         D(I) = XY(I)
    2 CONTINUE
      J = NN1
      DO 3 I = JJ , NPM
         J = J + 1
         XY(J) = D(I)
    3 CONTINUE
      DO 4 I = NN2 , NN
         J = J + 1
         XY(J) = D(I)
    4 CONTINUE
    5 NEQS = NSCON + 2*NPCON
      IF ( NEQS.EQ.0 ) GOTO 7
      NNEQ = N1C + N2C
      NEQS = NNEQ - NEQS + 1
C     Compute INV(ZA)E1.
      DO 6 I = NEQS , NNEQ
         XY(I) = (0.0D0,0.0D0)
    6 CONTINUE
    7 CALL SOLVES(ZA,IIP,XY,N1C,1,NNP,NN1,MMP,MM1,13,IFL,LD,D,IFAIL,
     &            DEBUG)
      IF(IFAIL.NE.0) RETURN
      NI = 0
      NPB = NPBL
C     Compute E2-ZC(INV(ZA)E1).
      DO 10 JJ = 1 , NBBL
         IF ( JJ.EQ.NBBL ) NPB = NLBL
         IF ( ICASX.GT.1 ) READ (CHTMP5) ((ZC(I,J),I=1,N1C),J=1,NPB)
         II = N1C + NI
         DO 9 I = 1 , NPB
            SUM = (0.0D0,0.0D0)
            DO 8 J = 1 , N1C
               SUM = SUM + ZC(J,I)*XY(J)
    8       CONTINUE
            J = II + I
            XY(J) = XY(J) - SUM
    9    CONTINUE
         NI = NI + NPBL
   10 CONTINUE
C     Following line from NEC81.
C     OPEN (CHTMP5,FORM='UNFORMATTED')
      REWIND CHTMP5
      JJ = N1C + 1
C     Compute INV(ZD)(E2-ZC(INV(ZA)E1)) = I2.
      IF ( ICASX.GT.1 ) GOTO 11
      CALL SOLVE(N2C,ZD,IIP(JJ),XY(JJ),N2C,LD,D,DEBUG)
      GOTO 13
   11 IF ( ICASX.EQ.4 ) GOTO 12
      NI = N2C*N2C
      READ (CHTMP1) (ZB(J,1),J=1,NI)
      REWIND CHTMP1
      CALL SOLVE(N2C,ZB,IIP(JJ),XY(JJ),N2C,LD,D,DEBUG)
      GOTO 13
   12 NBLSYS = NBLSYM
      NPSYS = NPSYM
      NLSYS = NLSYM
      ICASS = ICASE
      NBLSYM = NBBL
      NPSYM = NPBL
      NLSYM = NLBL
      ICASE = 3
      REWIND CHTMP1
      REWIND CHTMP6
      CALL LTSOLV(ZB,N2C,IIP(JJ),XY(JJ),N2C,1,11,16,LD,D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      REWIND CHTMP1
C     Following line from NEC81.
C     OPEN (CHTMP6,FORM='UNFORMATTED')
      REWIND CHTMP6
      NBLSYM = NBLSYS
      NPSYM = NPSYS
      NLSYM = NLSYS
      ICASE = ICASS
   13 NI = 0
      NPB = NPBL
C     Compute INV(ZA)E1-(INV(ZA)ZB)I2 = I1.
      DO 16 JJ = 1 , NBBL
         IF ( JJ.EQ.NBBL ) NPB = NLBL
         IF ( ICASX.GT.1 ) READ (CHTMP4) ((ZB(I,J),I=1,N1C),J=1,NPB)
         II = N1C + NI
         DO 15 I = 1 , N1C
            SUM = (0.0D0,0.0D0)
            DO 14 J = 1 , NPB
               JP = II + J
               SUM = SUM + ZB(I,J)*XY(JP)
   14       CONTINUE
            XY(I) = XY(I) - SUM
   15    CONTINUE
         NI = NI + NPBL
   16 CONTINUE
C     Following line from NEC81.
C     OPEN (CHTMP4,FORM='UNFORMATTED')
      REWIND CHTMP4
      IF ( NN1.EQ.NN .OR. MM1.EQ.0 ) GOTO 20
C     Reorder current array.
      DO 17 I = NN2 , NPM
         D(I) = XY(I)
   17 CONTINUE
      JJ = N1C + 1
      J = NN1
      DO 18 I = JJ , NPM
         J = J + 1
         XY(J) = D(I)
   18 CONTINUE
      DO 19 I = NN2 , N1C
         J = J + 1
         XY(J) = D(I)
   19 CONTINUE
   20 IF ( NSCON.EQ.0 ) GOTO 22
      J = NEQS - 1
      DO 21 I = 1 , NSCON
         J = J + 1
         JJ = ISCON(I)
         XY(JJ) = XY(J)
   21 CONTINUE
 
   22 RETURN
 
      END
