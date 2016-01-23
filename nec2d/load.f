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
 
      SUBROUTINE LOAD( LDTYP , LDTAG , LDTAGF , LDTAGT , ZLR , ZLI ,
     &                 ZLC , LD , NLOAD , ZARRAY , N1 , N2 , N , 
     &                 NP , M1 , ITAG , WLAM , SI , BI , IFAIL ,
     &                 DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* LOAD calculates the impedance of specified segments for various    **
C* types of loading.                                                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - LDTYP                                                     **
C* INPUT  - LDTAG                                                     **
C* INPUT  - LDTAGF                                                    **
C* INPUT  - LDTAGT                                                    **
C* INPUT  - ZLR                                                       **
C* INPUT  - ZLI                                                       **
C* INPUT  - ZLC                                                       **
C* INPUT  - LD                                                        **
C* INPUT  - NLOAD                                                     **
C* OUTPUT - ZARRAY                                                    **
C* INPUT  - N1                                                        **
C* INPUT  - N2                                                        **
C* INPUT  - N                                                         **
C* INPUT  - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - ITAG                                                      **
C* INPUT  - WLAM                                                      **
C* INPUT  - SI                                                        **
C* PASSED - BI                                                        **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       PRNT    ZINT                                           **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      COMPLEX*16 ZINT
      EXTERNAL PRNT , ZINT

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER LD , NLOAD , N1 , N2 , N , NP , M1 , IFAIL , DEBUG
      INTEGER LDTYP(1) , LDTAG(1) , LDTAGF(1) , LDTAGT(1) , 
     &        ITAG(2*LD)
      REAL*8 WLAM 
      REAL*8  ZLR(1) , ZLI(1) , ZLC(1) , SI(LD) , BI(LD)
      COMPLEX*16 ZARRAY(LD)

C     Local variables.
      INTEGER I , ICHK , ISTEP , IWARN , JUMP , L1 , L2 , LDTAGS , NOP
      REAL*8 RMU
      REAL*8 TPCJX(2)
      COMPLEX*16 ZT , TPCJ

C     Equivalent storage.
      EQUIVALENCE (TPCJ,TPCJX)
 
C     Data initialisation.
      DATA TPCJX /0.0D0 , 1.883651567308853277340992D9/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING LOAD'
      ENDIF

      WRITE (CHRSLT,26)

C     Initialize D array, used for temporary storage of loading
C     information.

      DO 1 I = N2 , N
         ZARRAY(I) = DCMPLX(0.0D0,0.0D0)
    1 CONTINUE
      IWARN = 0

C     Cycle over loading cards

      ISTEP = 0
    2 ISTEP = ISTEP + 1
      IF ( ISTEP.LE.NLOAD ) GOTO 6
      IF ( IWARN.EQ.1 ) WRITE (CHRSLT,27)
      IF ( N1+2*M1.GT.0 ) GOTO 5
      NOP = N/NP
      IF ( NOP.EQ.1 ) GOTO 5
      DO 4 I = 1 , NP
         ZT = ZARRAY(I)
         L1 = I
         DO 3 L2 = 2 , NOP
            L1 = L1 + NP
            ZARRAY(L1) = ZT
    3    CONTINUE
    4 CONTINUE
    5 RETURN
    6 IF ( LDTYP(ISTEP).LE.5 ) GOTO 7
      WRITE (CHRSLT,28) LDTYP(ISTEP)
      IFAIL=24
      RETURN
    7 LDTAGS = LDTAG(ISTEP)
      JUMP = LDTYP(ISTEP) + 1
      ICHK = 0

C     Search segments for proper ITAGS.

      L1 = N2
      L2 = N
      IF ( LDTAGS.NE.0 ) GOTO 8
      IF ( LDTAGF(ISTEP).EQ.0 .AND. LDTAGT(ISTEP).EQ.0 ) GOTO 8
      L1 = LDTAGF(ISTEP)
      L2 = LDTAGT(ISTEP)
      IF ( L1.GT.N1 ) GOTO 8
      WRITE (CHRSLT,30)
      IFAIL=25
      RETURN
    8 DO 18 I = L1 , L2
         IF ( LDTAGS.EQ.0 ) GOTO 9
         IF ( LDTAGS.NE.ITAG(I) ) GOTO 18
         IF ( LDTAGF(ISTEP).EQ.0 ) GOTO 9
         ICHK = ICHK + 1
         IF ( ICHK.GE.LDTAGF(ISTEP) .AND. ICHK.LE.LDTAGT(ISTEP) )
     &        GOTO 10
         GOTO 18
    9    ICHK = 1

C     Calculation of LAMDA*IMPED. per unit length, jump to appropriate
C     section for loading type.

   10    GOTO (11,12,13,14,15,16) , JUMP
   11    ZT = ZLR(ISTEP)/SI(I) + TPCJ*ZLI(ISTEP)/(SI(I)*WLAM)
         IF ( DABS(ZLC(ISTEP)).GT.1.D-20 ) ZT = ZT + 
     &        WLAM/(TPCJ*SI(I)*ZLC(ISTEP))
         GOTO 17
   12    ZT = TPCJ*SI(I)*ZLC(ISTEP)/WLAM
         IF ( DABS(ZLI(ISTEP)).GT.1.D-20 ) ZT = ZT + SI(I)
     &        *WLAM/(TPCJ*ZLI(ISTEP))
         IF ( DABS(ZLR(ISTEP)).GT.1.D-20 ) ZT = ZT + SI(I)/ZLR(ISTEP)
         ZT = 1.0D0/ZT
         GOTO 17
   13    ZT = ZLR(ISTEP)*WLAM + TPCJ*ZLI(ISTEP)
         IF ( DABS(ZLC(ISTEP)).GT.1.0D-20 ) ZT = ZT + 
     &        1.0D0/(TPCJ*SI(I)*SI(I)*ZLC(ISTEP))
         GOTO 17
   14    ZT = TPCJ*SI(I)*SI(I)*ZLC(ISTEP)
         IF ( DABS(ZLI(ISTEP)).GT.1.0D-20 ) 
     &        ZT = ZT + 1.0D0/(TPCJ*ZLI(ISTEP))
         IF ( DABS(ZLR(ISTEP)).GT.1.0D-20 ) 
     &        ZT = ZT + 1.0D0/(ZLR(ISTEP)*WLAM)
         ZT = 1.0D0/ZT
         GOTO 17
   15    ZT = DCMPLX(ZLR(ISTEP),ZLI(ISTEP))/SI(I)
         GOTO 17
C     IDF - fix from NEC81 for wire permeability.
C15   ZT=ZINT(ZLR(ISTEP)*WLAM,BI(I))
   16    RMU = ZLI(ISTEP)
         IF ( RMU.EQ.0.0D0 ) RMU = 1.0D0
         ZT = ZINT(ZLR(ISTEP)*WLAM,RMU,BI(I),DEBUG)
C     End of fix.
   17    IF ( (DABS(DBLE(ZARRAY(I)))+DABS(DIMAG(ZARRAY(I))))
     &        .GT.1.D-20 ) IWARN = 1
 
 
         ZARRAY(I) = ZARRAY(I) + ZT
   18 CONTINUE
      IF ( ICHK.NE.0 ) GOTO 19
      WRITE (CHRSLT,29) LDTAGS
      IFAIL=26
      RETURN

C     Printing the segment loading data, jump to proper print.

   19 GOTO (20,21,22,23,24,25) , JUMP
   20 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),ZLR(ISTEP),ZLI(ISTEP)
     &          ,ZLC(ISTEP),0.0D0,0.0D0,0.0D0,' SERIES ',DEBUG)
      GOTO 2
   21 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),ZLR(ISTEP),ZLI(ISTEP)
     &          ,ZLC(ISTEP),0.0D0,0.0D0,0.0D0,'PARALLEL',DEBUG)
      GOTO 2
   22 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),ZLR(ISTEP),ZLI(ISTEP)
     &          ,ZLC(ISTEP),0.0D0,0.0D0,0.0D0,' SERIES (PER METER) ',
     &          DEBUG)
      GOTO 2
   23 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),ZLR(ISTEP),ZLI(ISTEP)
     &          ,ZLC(ISTEP),0.0D0,0.0D0,0.0D0,'PARALLEL (PER METER)',
     &          DEBUG)
      GOTO 2
   24 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),0.0D0,0.0D0,0.0D0,
     &          ZLR(ISTEP),ZLI(ISTEP),0.0D0,'FIXED IMPEDANCE ',DEBUG)
      GOTO 2
   25 CALL PRNT(LDTAGS,LDTAGF(ISTEP),LDTAGT(ISTEP),0.0D0,0.0D0,0.0D0,
     &          0.0D0,0.0D0,ZLR(ISTEP),'  WIRE  ',DEBUG)
      GOTO 2
 
   26 FORMAT (//,7X,'LOCATION',10X,'RESISTANCE',3X,'INDUCTANCE',2X,
     &        'CAPACITANCE',7X,'IMPEDANCE (OHMS)',5X,'CONDUCTIVITY',4X,
     &        'TYPE',/,4X,'ITAG',' FROM THRU',10X,'OHMS',8X,'HENRYS',7X,
     &        'FARADS',8X,'REAL',6X,'IMAGINARY',4X,'MHOS/METER')
   27 FORMAT (/,10X,
     &'NOTE, SOME OF THE ABOVE SEGMENTS HAVE BEEN LOADED TWICE - IMPEDAN
     &CES ADDED')
   28 FORMAT (/,10X,'IMPROPER LOAD TYPE CHOOSEN, REQUESTED TYPE IS ',I3)
   29 FORMAT (/,10X,'LOADING DATA CARD ERROR, NO SEGMENT HAS AN ITAG = '
     &        ,I5)
   30 FORMAT (
     & ' ERROR - LOADING MAY NOT BE ADDED TO SEGMENTS IN N.G.F. SECTION'
     & )
 
      END
