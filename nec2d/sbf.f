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

      SUBROUTINE SBF( I , IS , AA , BB , CC , LD , ICON1 , ICON2 , 
     &                SI , BI , JMAX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute component of basis function I on segment IS.               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - I                                                         **
C* INPUT  - IS                                                        **
C* OUTPUT - AA                                                        **
C* OUTPUT - BB                                                        **
C* OUTPUT - CC                                                        **
C* INPUT  - LD                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - JMAX                                                      **
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
C* called by   TRIO                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter defintions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER I , IS , LD , JMAX , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD)
      REAL*8 AA , BB , CC 
      REAL*8 SI(LD) , BI(LD)

C     Local variables.
      INTEGER IEND , JCOX , JEND , JJSNO , JUNE , NJUN1 , NJUN2
      REAL*8 AJ , AP , CD , CDH , DD , OMC , PI , PM , PP , QM , QP , 
     &       SD , SDH , SSIG , XXI , EULER , OO24 , OO720
 
C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
      DATA EULER /0.577215664901532860606512D0/
      DATA OO24 /0.04166666666666666666666666D0/
      DATA OO720 /0.001388888888888888888888888D0/

 
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SBF'
      ENDIF
      
      AA = 0.0D0
      BB = 0.0D0
      CC = 0.0D0
      JUNE = 0
      JJSNO = 0
      PP = 0.0D0
      JCOX = ICON1(I)
      IF ( JCOX.GT.PATOFF ) JCOX = I
      JEND = -1
      IEND = -1
      SSIG = -1.0D0
      IF ( JCOX ) 1 , 11 , 2
    1 JCOX = -JCOX
      GOTO 3
    2 SSIG = -SSIG
      JEND = -JEND
    3 JJSNO = JJSNO + 1
      IF ( JJSNO.GE.JMAX ) GOTO 24
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
      IF ( JCOX.NE.IS ) GOTO 6
      AA = AJ/SD*SSIG
      BB = AJ/(2.0D0*CDH)
      CC = -AJ/(2.0D0*SDH)*SSIG
      JUNE = IEND
    6 IF ( JCOX.EQ.I ) GOTO 9
      IF ( JEND.EQ.1 ) GOTO 7
      JCOX = ICON1(JCOX)
      GOTO 8
    7 JCOX = ICON2(JCOX)
    8 IF ( IABS(JCOX).EQ.I ) GOTO 10
      IF ( JCOX ) 1 , 24 , 2
    9 IF ( JCOX.EQ.IS ) BB = -BB
   10 IF ( IEND.EQ.1 ) GOTO 12
   11 PM = -PP
      PP = 0.0D0
      NJUN1 = JJSNO
      JCOX = ICON2(I)
      IF ( JCOX.GT.PATOFF ) JCOX = I
      JEND = 1
      IEND = 1
      SSIG = -1.0D0
      IF ( JCOX ) 1 , 12 , 2
   12 NJUN2 = JJSNO - NJUN1
      DD = PI*SI(I)
      SDH = DSIN(DD)
      CDH = DCOS(DD)
      SD = 2.0D0*SDH*CDH
      CD = CDH*CDH - SDH*SDH
      IF ( DD.GT.0.015D0 ) GOTO 13
      OMC = 4.0D0*DD*DD
      OMC = ((OO720*OMC-OO24)*OMC+0.5D0)*OMC
      GOTO 14
   13 OMC = 1.0D0 - CD
   14 AP = 1.0D0/(DLOG(1.0D0/(PI*BI(I)))-EULER)
      AJ = AP
      IF ( NJUN1.EQ.0 ) GOTO 19
      IF ( NJUN2.EQ.0 ) GOTO 21
      QP = SD*(PM*PP+AJ*AP) + CD*(PM*AP-PP*AJ)
      QM = (AP*OMC-PP*SD)/QP
      QP = -(AJ*OMC+PM*SD)/QP
      IF ( JUNE ) 15 , 18 , 16
   15 AA = AA*QM
      BB = BB*QM
      CC = CC*QM
      GOTO 17
   16 AA = -AA*QP
      BB = BB*QP
      CC = -CC*QP
   17 IF ( I.NE.IS ) RETURN
   18 AA = AA - 1.0D0
      BB = BB + (AJ*QM+AP*QP)*SDH/SD
      CC = CC + (AJ*QM-AP*QP)*CDH/SD
      RETURN
   19 IF ( NJUN2.EQ.0 ) GOTO 23
      QP = PI*BI(I)
      XXI = QP*QP
      XXI = QP*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
      QP = -(OMC+XXI*SD)/(SD*(AP+XXI*PP)+CD*(XXI*AP-PP))
      IF ( JUNE.NE.1 ) GOTO 20
      AA = -AA*QP
      BB = BB*QP
      CC = -CC*QP
      IF ( I.NE.IS ) RETURN
   20 AA = AA - 1.0D0
      DD = CD - XXI*SD
      BB = BB + (SDH+AP*QP*(CDH-XXI*SDH))/DD
      CC = CC + (CDH+AP*QP*(SDH+XXI*CDH))/DD
      RETURN
   21 QM = PI*BI(I)
      XXI = QM*QM
      XXI = QM*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
      QM = (OMC+XXI*SD)/(SD*(AJ-XXI*PM)+CD*(PM+XXI*AJ))
      IF ( JUNE.NE.-1 ) GOTO 22
      AA = AA*QM
      BB = BB*QM
      CC = CC*QM
      IF ( I.NE.IS ) RETURN
   22 AA = AA - 1.0D0
      DD = CD - XXI*SD
      BB = BB + (AJ*QM*(CDH-XXI*SDH)-SDH)/DD
      CC = CC + (CDH-AJ*QM*(SDH+XXI*CDH))/DD
      RETURN
   23 AA = -1.0D0
      QP = PI*BI(I)
      XXI = QP*QP
      XXI = QP*(1.0D0-0.5D0*XXI)/(1.0D0-XXI)
      CC = 1.0D0/(CDH-XXI*SDH)
      RETURN
   24 WRITE (CHRSLT,25) I
 
      IFAIL=33
      RETURN
 
   25 FORMAT (' SBF - SEGMENT CONNECTION ERROR FOR SEGMENT',I5)
 
      END
