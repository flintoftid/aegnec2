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

      SUBROUTINE TBF( I , ICAP , LD , ICON1 , ICON2 , SI , BI , 
     &                JMAX , JSNO , JCO , AX , BX , CX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C*    Compute the basis function for segment I.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - I                                                         **
C* INPUT  - ICAP                                                      **
C* INPUT  - LD                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - JMAX                                                      **
C* OUTPUT - JSNO                                                      **
C* OUTPUT - JCO                                                       **
C* OUTPUT - AX                                                        **
C* OUTPUT - BX                                                        **
C* OUTPUT - CX                                                        **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CABC    CMNGF   QDSRC                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER LD , I , ICAP , JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , JCO(JMAX)
      REAL*8 SI(LD) , BI(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX) 

C     Local variables.
      INTEGER IEND , JCOX , JEND , JSNOP , NJUN1 , NJUN2
      REAL*8 AJ , AP , CD , CDH , DD , OMC , PI , PM , PP , QM , 
     &       QP , SD , SDH , SSIG , XXI , EULER , OO24 , OO720
 
C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA EULER /0.577215664901532860606512D0/
      DATA OO24 /0.04166666666666666666666667D0/
      DATA OO720 /0.001388888888888888888888889D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING TBF'
      ENDIF

      JSNO = 0
      PP = 0.0D0
      JCOX = ICON1(I)
      IF ( JCOX.GT.PATOFF ) JCOX = I
      JEND = -1
      IEND = -1
      SSIG = -1.0D0
      IF ( JCOX ) 1 , 10 , 2
    1 JCOX = -JCOX
      GOTO 3
    2 SSIG = -SSIG
      JEND = -JEND
    3 JSNO = JSNO + 1
      IF ( JSNO.GE.JMAX ) GOTO 28
      JCO(JSNO) = JCOX
      DD = PI*SI(JCOX)
      SDH = DSIN(DD)
      CDH = DCOS(DD)
      SD = 2.0D0*SDH*CDH
      IF ( DD.GT.0.015D0 ) GOTO 4
      OMC = 4.0D0*DD*DD
      OMC = ((OO720*OMC-OO24)*OMC+0.5D0)*OMC
      GOTO 5
    4 OMC = 1.0D0 - CDH*CDH + SDH*SDH
    5 AJ = 1.0D0/(DLOG(1.0D0/(PI*BI(JCOX)))-EULER)
      PP = PP - OMC/SD*AJ
      AX(JSNO) = AJ/SD*SSIG
      BX(JSNO) = AJ/(2.0D0*CDH)
      CX(JSNO) = -AJ/(2.0D0*SDH)*SSIG
      IF ( JCOX.EQ.I ) GOTO 8
      IF ( JEND.EQ.1 ) GOTO 6
      JCOX = ICON1(JCOX)
      GOTO 7
    6 JCOX = ICON2(JCOX)
    7 IF ( IABS(JCOX).EQ.I ) GOTO 9
      IF ( JCOX ) 1 , 28 , 2
    8 BX(JSNO) = -BX(JSNO)
    9 IF ( IEND.EQ.1 ) GOTO 11
   10 PM = -PP
      PP = 0.0D0
      NJUN1 = JSNO
      JCOX = ICON2(I)
      IF ( JCOX.GT.PATOFF ) JCOX = I
      JEND = 1
      IEND = 1
      SSIG = -1.0D0
      IF ( JCOX ) 1 , 11 , 2
   11 NJUN2 = JSNO - NJUN1
      JSNOP = JSNO + 1
      JCO(JSNOP) = I
      DD = PI*SI(I)
      SDH = DSIN(DD)
      CDH = DCOS(DD)
      SD = 2.0D0*SDH*CDH
      CD = CDH*CDH - SDH*SDH
      IF ( DD.GT.0.015D0 ) GOTO 12
      OMC = 4.0D0*DD*DD
      OMC = ((OO720*OMC-OO24)*OMC+0.5D0)*OMC
      GOTO 13
   12 OMC = 1.0D0 - CD
   13 AP = 1.0D0/(DLOG(1.0D0/(PI*BI(I)))-EULER)
      AJ = AP
      IF ( NJUN1.EQ.0 ) GOTO 16
      IF ( NJUN2.EQ.0 ) GOTO 20
      QP = SD*(PM*PP+AJ*AP) + CD*(PM*AP-PP*AJ)
      QM = (AP*OMC-PP*SD)/QP
      QP = -(AJ*OMC+PM*SD)/QP
      BX(JSNOP) = (AJ*QM+AP*QP)*SDH/SD
      CX(JSNOP) = (AJ*QM-AP*QP)*CDH/SD
      DO 14 IEND = 1 , NJUN1
         AX(IEND) = AX(IEND)*QM
         BX(IEND) = BX(IEND)*QM
         CX(IEND) = CX(IEND)*QM
   14 CONTINUE
      JEND = NJUN1 + 1
      DO 15 IEND = JEND , JSNO
         AX(IEND) = -AX(IEND)*QP
         BX(IEND) = BX(IEND)*QP
         CX(IEND) = -CX(IEND)*QP
   15 CONTINUE
      GOTO 27
   16 IF ( NJUN2.EQ.0 ) GOTO 24
      IF ( ICAP.NE.0 ) GOTO 17
      XXI = 0.0D0
      GOTO 18
   17 QP = PI*BI(I)
      XXI = QP*QP
      XXI = QP*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
   18 QP = -(OMC+XXI*SD)/(SD*(AP+XXI*PP)+CD*(XXI*AP-PP))
      DD = CD - XXI*SD
      BX(JSNOP) = (SDH+AP*QP*(CDH-XXI*SDH))/DD
      CX(JSNOP) = (CDH+AP*QP*(SDH+XXI*CDH))/DD
      DO 19 IEND = 1 , NJUN2
         AX(IEND) = -AX(IEND)*QP
         BX(IEND) = BX(IEND)*QP
         CX(IEND) = -CX(IEND)*QP
   19 CONTINUE
      GOTO 27
   20 IF ( ICAP.NE.0 ) GOTO 21
      XXI = 0.0D0
      GOTO 22
   21 QM = PI*BI(I)
      XXI = QM*QM
      XXI = QM*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
   22 QM = (OMC+XXI*SD)/(SD*(AJ-XXI*PM)+CD*(PM+XXI*AJ))
      DD = CD - XXI*SD
      BX(JSNOP) = (AJ*QM*(CDH-XXI*SDH)-SDH)/DD
      CX(JSNOP) = (CDH-AJ*QM*(SDH+XXI*CDH))/DD
      DO 23 IEND = 1 , NJUN1
         AX(IEND) = AX(IEND)*QM
         BX(IEND) = BX(IEND)*QM
         CX(IEND) = CX(IEND)*QM
   23 CONTINUE
      GOTO 27
   24 BX(JSNOP) = 0.0D0
      IF ( ICAP.NE.0 ) GOTO 25
      XXI = 0.0D0
      GOTO 26
   25 QP = PI*BI(I)
      XXI = QP*QP
      XXI = QP*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
   26 CX(JSNOP) = 1.0D0/(CDH-XXI*SDH)
   27 JSNO = JSNOP
      AX(JSNO) = -1.0D0
      RETURN
   28 WRITE (CHRSLT,29) I
 
      IFAIL=33
      RETURN
 
   29 FORMAT (' TBF - SEGMENT CONNECTION ERROR FOR SEGMENT',I5)
 
      END
