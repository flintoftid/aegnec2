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

      SUBROUTINE GFIL( IPRT , LD , SALP , IRESRV , CM , NLODF , 
     &                 ZARRAY , KCOM , COM  , IP , EPSR , SIG , SCRWLT , 
     &                 SCRWRT , FMHZ , IPSYM , N1 , N2 , N , NP , M1 , 
     &                 M2 , M, MP , ICON1 , ICON2 , ITAG , WLAM , X , 
     &                 Y , Z , SI , BI , ALP , BET , T2X , T2Y , T2Z , 
     &                 JOBNAM , FILBUF , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* GFIL reads the NGF file.                                           **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - IPRT                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
C* INPUT  - IRESRV                                                    **
C* OUTPUT - CM                                                        **
C* OUTPUT - NLODF                                                     **
C* OUTPUT - ZARRAY                                                    **
C* OUTPUT - KCOM                                                      **
C* OUTPUT - COM                                                       **
C* OUTPUT - IP                                                        **
C* OUTPUT - EPSR                                                      **
C* OUTPUT - SIG                                                       **
C* OUTPUT - SCRWLT                                                    **
C* OUTPUT - SCRWRT                                                    **
C* OUTPUT - FMHZ                                                      **
C* OUTPUT - IPSYM                                                     **
C* OUTPUT - N1                                                        **
C* OUTPUT - N2                                                        **
C* OUTPUT - N                                                         **
C* OUTPUT - NP                                                        **
C* OUTPUT - M1                                                        **
C* OUTPUT - M2                                                        **
C* OUTPUT - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - ICON1                                                     **
C* OUTPUT - ICON2                                                     **
C* OUTPUT - ITAG                                                      **
C* OUTPUT - WLAM                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - SI                                                        **
C* OUTPUT - BI                                                        **
C* OUTPUT - ALP                                                       **
C* OUTPUT - BET                                                       **
C* OUTPUT - T2X                                                       **
C* OUTPUT - T2Y                                                       **
C* OUTPUT - T2Z                                                       **
C* PASSED - JOBNAM                                                    **
C* INPUT  - FILBUF                                                    **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    AR1     AR2     AR3     DXA     DYA     EPSCF   ICASE  **
C*             IMAT    IPERF   KSYMP   NBLOKS  NBLSYM  NLAST   NLSYM  **
C*             NPBLK   NPSYM   NRADL   NXA     NYA     SSX     XSA    **
C*             YSA                                                    **
C* uses value  ICASE   IMAT    IPERF   KSYMP   NBLSYM  NPSYM          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN  BLCKOT  FILSUF                                 **
C* called by   DATAGN                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines
      EXTERNAL BLCKIN , BLCKOT , FILSUF

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF
      CHARACTER*78 COM(5)
      INTEGER IPRT , LD , IRESRV , NLODF , KCOM , IPSYM , N1 , N2 , N , 
     &        NP , M1 , M2 , M , MP , IFAIL , DEBUG
      INTEGER IP(2*LD) , ICON1(2*LD) , ICON2(2*LD) , 
     &        ITAG(2*LD)
      REAL*8 EPSR , SIG , SCRWLT , SCRWRT , FMHZ , WLAM 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , BI(LD) , 
     &       ALP(LD) , BET(LD) , T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 CM(IRESRV) , ZARRAY(LD)

C     Local variables.
      INTEGER I , IOP , IOUT , J , K , NBL2 , NNEQ , NOP , NNPEQ 
      REAL*8 DX , XI , YI , ZI

C     Common storage.
      INCLUDE 'gnd.inc'
      INCLUDE 'ggrid.inc'
      INCLUDE 'matpar.inc'
      INCLUDE 'smat.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING GFIL'
      ENDIF
 
      CALL FILSUF(JOBNAM,FILBUF,0,SUFNGF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      OPEN (UNIT=CHNGFL,FILE=FILBUF,FORM='UNFORMATTED',
     &      STATUS='OLD')
      REWIND CHNGFL
      READ (CHNGFL) N1 , NP , M1 , MP , WLAM , FMHZ , IPSYM , KSYMP , 
     &            IPERF , NRADL , EPSR , SIG , SCRWLT , SCRWRT , NLODF , 
     &            KCOM
      N = N1
      M = M1
      N2 = N1 + 1
      M2 = M1 + 1
      IF ( N1.EQ.0 ) GOTO 2
C     Read seg. data and convert back to end coord. in units of meters.
      READ (CHNGFL) (X(I),I=1,N1) , (Y(I),I=1,N1) , (Z(I),I=1,N1)
      READ (CHNGFL) (SI(I),I=1,N1) , (BI(I),I=1,N1) , (ALP(I),I=1,N1)
      READ (CHNGFL) (BET(I),I=1,N1) , (SALP(I),I=1,N1)
      READ (CHNGFL) (ICON1(I),I=1,N1) , (ICON2(I),I=1,N1)
      READ (CHNGFL) (ITAG(I),I=1,N1)
      IF ( NLODF.NE.0 ) READ (CHNGFL) (ZARRAY(I),I=1,N1)
      DO 1 I = 1 , N1
         XI = X(I)*WLAM
         YI = Y(I)*WLAM
         ZI = Z(I)*WLAM
         DX = SI(I)*0.5D0*WLAM
         X(I) = XI - ALP(I)*DX
         Y(I) = YI - BET(I)*DX
         Z(I) = ZI - SALP(I)*DX
         SI(I) = XI + ALP(I)*DX
         ALP(I) = YI + BET(I)*DX
         BET(I) = ZI + SALP(I)*DX
         BI(I) = BI(I)*WLAM
    1 CONTINUE
    2 IF ( M1.EQ.0 ) GOTO 4
      J = LD - M1 + 1
C     Read patch data and convert to meters.
      READ (CHNGFL) (X(I),I=J,LD) , (Y(I),I=J,LD) , (Z(I),I=J,LD)
      READ (CHNGFL) (SI(I),I=J,LD) , (BI(I),I=J,LD) , (ALP(I),I=J,LD)
      READ (CHNGFL) (BET(I),I=J,LD) , (SALP(I),I=J,LD)
C     Error corrected 11/20/89.
      READ (CHNGFL) (T2X(I),I=J,LD) , (T2Y(I),I=J,LD)
      READ (CHNGFL) (T2Z(I),I=J,LD)
C     READ (CHNGFL) (ICON1(I),I=J,LD),(ICON2(I),I=J,LD)
C     READ (CHNGFL) (ITAG(I),I=J,LD)

      DX = WLAM*WLAM
      DO 3 I = J , LD
         X(I) = X(I)*WLAM
         Y(I) = Y(I)*WLAM
         Z(I) = Z(I)*WLAM
         BI(I) = BI(I)*DX
    3 CONTINUE
    4 READ (CHNGFL) ICASE , NBLOKS , NPBLK , NLAST , NBLSYM , NPSYM , 
     &            NLSYM , IMAT
      IF ( IPERF.EQ.2 ) READ (CHNGFL) AR1 , AR2 , AR3 , EPSCF , DXA , 
     &                              DYA , XSA , YSA , NXA , NYA
      NNEQ = N1 + 2*M1
      NNPEQ = NP + 2*MP
      NOP = NNEQ/NNPEQ
      IF ( NOP.GT.1 ) READ (CHNGFL) ((SSX(I,J),I=1,NOP),J=1,NOP)
      READ (CHNGFL) (IP(I),I=1,NNEQ) , COM
C     Read matrix A and write tape 13 for out of core.
      IF ( ICASE.GT.2 ) GOTO 5
      IOUT = NNEQ*NNPEQ
      READ (CHNGFL) (CM(I),I=1,IOUT)
      GOTO 11
    5 REWIND CHTMP3
      IF ( ICASE.NE.4 ) GOTO 7
      IOUT = NNPEQ*NNPEQ
      DO 6 K = 1 , NOP
         READ (CHNGFL) (CM(J),J=1,IOUT)
         WRITE (CHTMP3) (CM(J),J=1,IOUT)
    6 CONTINUE
      GOTO 10
    7 IOUT = NPSYM*NNPEQ*2
      NBL2 = 2*NBLSYM
      DO 9 IOP = 1 , NOP
         DO 8 I = 1 , NBL2
            CALL BLCKIN(CM,CHNGFL,1,IOUT,1,206,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            CALL BLCKOT(CM,13,1,IOUT,DEBUG)
    8    CONTINUE
    9 CONTINUE
   10 REWIND CHTMP3
   11 REWIND CHNGFL

      WRITE (CHRSLT,17)
      WRITE (CHRSLT,15)
      WRITE (CHRSLT,15)
      WRITE (CHRSLT,18)
      WRITE (CHRSLT,19) N1 , M1
      IF ( NOP.GT.1 ) WRITE (CHRSLT,20) NOP
      WRITE (CHRSLT,21) IMAT , ICASE
      IF ( ICASE.LT.3 ) GOTO 12
      NBL2 = NNEQ*NNPEQ
      WRITE (CHRSLT,22) NBL2
   12 WRITE (CHRSLT,23) FMHZ
      IF ( KSYMP.EQ.2 .AND. IPERF.EQ.1 ) WRITE (CHRSLT,24)
      IF ( KSYMP.EQ.2 .AND. IPERF.EQ.0 ) WRITE (CHRSLT,28)
      IF ( KSYMP.EQ.2 .AND. IPERF.EQ.2 ) WRITE (CHRSLT,29)
      IF ( KSYMP.EQ.2 .AND. IPERF.NE.1 ) WRITE (CHRSLT,25) EPSR , SIG
      WRITE (CHRSLT,18)
      DO 13 J = 1 , KCOM
         WRITE (CHRSLT,16) COM(J)
   13 CONTINUE
      WRITE (CHRSLT,18)
      WRITE (CHRSLT,15)
      WRITE (CHRSLT,15)
      WRITE (CHRSLT,17)
C     From NEC81 - close file before return.
      IF ( IPRT.EQ.0 ) THEN
         CLOSE (CHNGFL)
         RETURN
      ENDIF
C
      WRITE (CHRSLT,26)
      DO 14 I = 1 , N1
         WRITE (CHRSLT,27) I , X(I) , Y(I) , Z(I) , SI(I) , ALP(I) , 
     &                     BET(I)
   14 CONTINUE
 
      RETURN
 
   15 FORMAT (5X,'**************************************************',
     &        '**********************************')
   16 FORMAT (5X,'** ',A,' **')
   17 FORMAT (////)
   18 FORMAT (5X,'**',80X,'**')
   19 FORMAT (5X,'** NUMERICAL GREEN''S FUNCTION',53X,'**',/,5X,
     &        '** NO. SEGMENTS =',I4,10X,'NO. PATCHES =',I4,34X,'**')
   20 FORMAT (5X,'** NO. SYMMETRIC SECTIONS =',I4,51X,'**')
   21 FORMAT (5X,'** N.G.F. MATRIX -  CORE STORAGE =',I7,
     &        ' COMPLEX NUMBERS,  CASE',I2,16X,'**')
   22 FORMAT (5X,'**',19X,'MATRIX SIZE =',I7,' COMPLEX NUMBERS',25X,
     &        '**')
   23 FORMAT (5X,'** FREQUENCY =',1P,E12.5,' MHZ.',51X,'**')
   24 FORMAT (5X,'** PERFECT GROUND',65X,'**')
   25 FORMAT (5X,'** GROUND PARAMETERS - DIELECTRIC CONSTANT =',1P,
     &        E12.5,26X,'**',/,5X,'**',21X,'CONDUCTIVITY =',E12.5,
     &        ' MHOS/M.',25X,'**')
   26 FORMAT (39X,'NUMERICAL GREEN''S FUNCTION DATA',/,41X,
     &        'COORDINATES OF SEGMENT ENDS',/,51X,'(METERS)',/,5X,
     &        'SEG.',11X,'- - - END ONE - - -',26X,
     &        '- - - END TWO - - -',/,6X,'NO.',6X,'X',14X,'Y',14X,'Z',
     &        14X,'X',14X,'Y',14X,'Z')
   27 FORMAT (1X,I7,1P,6E15.6)
   28 FORMAT (5X,
     &        '** FINITE GROUND.  REFLECTION COEFFICIENT APPROXIMATION',
     &        27X,'**')
   29 FORMAT (5X,'** FINITE GROUND.  SOMMERFELD SOLUTION',44X,'**')
 
      END
