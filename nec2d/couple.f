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

      SUBROUTINE COUPLE( CCUR , LD , N , ITAG , WLAM , NSMAX , 
     &                   NVQD , NSANT , ISANT , VSANT , MXCOUP , 
     &                   MXCOU2 , NCOUP , ICOUP , NCTAG , NCSEG , Y11A , 
     &                   Y12A , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* COUPLE computes the maximum coupling between pairs of segments.    **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - CCUR                                                      **
C* INPUT  - LD                                                        **
C* PASSED - N                                                         **
C* PASSED - ITAG                                                      **
C* INPUT  - WLAM                                                      **
C* INPUT  - NSMAX                                                     **
C* INPUT  - NVQD                                                      **
C* INPUT  - NSANT                                                     **
C* INPUT  - ISANT                                                     **
C* INPUT  - VSANT                                                     **
C* INPUT  - MXCOUP                                                    **
C* INPUT  - MXCOU2                                                    **
C* INPUT  - NCOUP                                                     **
C* OUTPUT - ICOUP                                                     **
C* INPUT  - NCTAG                                                     **
C* INPUT  - NCSEG                                                     **
C* OUTPUT - Y11A                                                      **
C* OUTPUT - Y12A                                                      **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       DB10    ISEGNO                                         **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      INTEGER ISEGNO
      REAL*8 DB10
      EXTERNAL DB10 , ISEGNO

C     Dummy arguments.
      INTEGER LD , N , NSMAX , NVQD , NSANT , MXCOUP , MXCOU2 , NCOUP , 
     &        ICOUP , IFAIL , DEBUG
      INTEGER ITAG(2*LD) , ISANT(NSMAX) , NCTAG(MXCOUP) , 
     &        NCSEG(MXCOUP)
      REAL*8 WLAM 
      COMPLEX*16 CCUR(1) , VSANT(NSMAX) , Y11A(MXCOU2) , Y12A(MXCOU2) 

C     Local variables.
      INTEGER I , ISG1 , ISG2 , ITS1 , ITS2 , ITT1 , ITT2 , 
     &        J , J1 , J2 , K , L1 , NPM1
      REAL*8 C , DBC , GMAX
      COMPLEX*16 Y11 , Y12 , Y22 , YL , YIN , ZL , ZIN , RHO


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING COUPLE'
      ENDIF
 
      IF ( NSANT.NE.1 .OR. NVQD.NE.0 ) RETURN
      J = ISEGNO(NCTAG(ICOUP+1),NCSEG(ICOUP+1),LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.EQ.0) RETURN
      IF ( J.NE.ISANT(1) ) RETURN
      ICOUP = ICOUP + 1
      ZIN = VSANT(1)
      Y11A(ICOUP) = CCUR(J)*WLAM/ZIN
      L1 = (ICOUP-1)*(NCOUP-1)
      DO 1 I = 1 , NCOUP
         IF ( I.EQ.ICOUP ) GOTO 1
         K = ISEGNO(NCTAG(I),NCSEG(I),LD,ITAG,N,IFAIL,DEBUG)
         IF(IFAIL.EQ.0) RETURN
         L1 = L1 + 1
         Y12A(L1) = CCUR(K)*WLAM/ZIN
    1 CONTINUE
      IF ( ICOUP.LT.NCOUP ) RETURN
      WRITE (3,7)
      NPM1 = NCOUP - 1
      DO 6 I = 1 , NPM1
         ITT1 = NCTAG(I)
         ITS1 = NCSEG(I)
         ISG1 = ISEGNO(ITT1,ITS1,LD,ITAG,N,IFAIL,DEBUG)
         IF(IFAIL.EQ.0) RETURN
         L1 = I + 1
         DO 5 J = L1 , NCOUP
            ITT2 = NCTAG(J)
            ITS2 = NCSEG(J)
            ISG2 = ISEGNO(ITT2,ITS2,LD,ITAG,N,IFAIL,DEBUG)
            IF(IFAIL.EQ.0) RETURN
            J1 = J + (I-1)*NPM1 - 1
            J2 = I + (J-1)*NPM1
            Y11 = Y11A(I)
            Y22 = Y11A(J)
            Y12 = (0.5D0,0.0D0)*(Y12A(J1)+Y12A(J2))
            YIN = Y12*Y12
            DBC = CDABS(YIN)
            C = DBC/(2.0D0*DBLE(Y11)*DBLE(Y22)-DBLE(YIN))
            IF ( C.LT.0.0D0 .OR. C.GT.1.0D0 ) GOTO 4
            IF ( C.LT.0.01D0 ) GOTO 2
            GMAX = (1.0D0-DSQRT(1.0D0-C*C))/C
            GOTO 3
    2       GMAX = 0.5D0*(C+0.25D0*C*C*C)
    3       RHO = GMAX*DCONJG(YIN)/DBC
            YL = ((1.0D0-RHO)/(1.0D0+RHO)+1.0D0)*DBLE(Y22) - Y22
            ZL = 1.0D0/YL
            YIN = Y11 - YIN/(Y22+YL)
            ZIN = 1.0D0/YIN
            DBC = DB10(GMAX)
            WRITE (3,8) ITT1 , ITS1 , ISG1 , ITT2 , ITS2 , ISG2 , DBC , 
     &                  ZL , ZIN
            GOTO 5
    4       WRITE (3,9) ITT1 , ITS1 , ISG1 , ITT2 , ITS2 , ISG2 , C
    5    CONTINUE
    6 CONTINUE
 
      RETURN
 
    7 FORMAT (///,36X,'- - - ISOLATION DATA - - -',//,6X,
     &        '- - COUPLING BETWEEN - -',8X,'MAXIMUM',15X,
     &        '- - - FOR MAXIMUM COUPLING - - -',/,12X,'SEG.',14X,
     &        'SEG.',3X,'COUPLING',4X,'LOAD IMPEDANCE (2ND SEG.)',7X,
     &        'INPUT IMPEDANCE',/,2X,'TAG/SEG.',3X,'NO.',4X,'TAG/SEG.',
     &        3X,'NO.',6X,'(DB)',8X,'REAL',9X,'IMAG.',9X,'REAL',9X,
     &        'IMAG.')
    8 FORMAT (2(1X,I4,1X,I4,1X,I5,2X),F9.3,2X,1P,2(2X,E12.5,1X,E12.5))
    9 FORMAT (2(1X,I4,1X,I4,1X,I5,2X),
     &        '**ERROR** COUPLING IS NOT BETWEEN 0 AND 1. (=',1P,E12.5,
     &        ')')
 
      END
