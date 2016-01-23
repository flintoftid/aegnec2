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

      SUBROUTINE GFOUT( LD , SALP , IRESRV , CM , NLOAD , NLODF , 
     &                  ZARRAY , KCOM , COM , IP, EPSR , SIG , SCRWLT , 
     &                  SCRWRT , FMHZ , IPSYM , N1 , N2 , N , NP , M1 , 
     &                  M2 , M , MP , ICON1 , ICON2 , ITAG , ICONX , 
     &                  WLAM , X , Y , Z , SI , BI , ALP , BET , T2X , 
     &                  T2Y , T2Z , JOBNAM , FILBUF , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Write NGF file.                                                    **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - IRESRV                                                    **
C* OUTPUT - CM                                                        **
C* INPUT  - NLOAD                                                     **
C* INPUT  - NLODF                                                     **
C* INPUT  - ZARRAY                                                    **
C* INPUT  - KCOM                                                      **
C* INPUT  - COM                                                       **
C* INPUT  - IP                                                        **
C* INPUT  - EPSR                                                      **
C* INPUT  - SIG                                                       **
C* INPUT  - SCRWLT                                                    **
C* INPUT  - SCRWRT                                                    **
C* INPUT  - FMHZ                                                      **
C* INPUT  - IPSYM                                                     **
C* INPUT  - N1                                                        **
C* INPUT  - N2                                                        **
C* INPUT  - N                                                         **
C* INPUT  - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M2                                                        **
C* INPUT  - M                                                         **
C* INPUT  - MP                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - ITAG                                                      **
C* INPUT  - ICONX                                                     **
C* INPUT  - WLAM                                                      **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - ALP                                                       **
C* INPUT  - BET                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* PASSED - JOBNAM                                                    **
C* INPUT  - FILBUF                                                    **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  AR1     AR2     AR3     CH      CL      DXA     DYA    **
C*             EPSCF   FRATI   ICASE   ICASX   IFAR    IMAT    IPERF  **
C*             KSYMP   NBBL    NBBX    NBLOKS  NBLSYM  NLAST   NLBL   **
C*             NLBX    NLSYM   NPBL    NPBLK   NPBX    NPSYM   NRADL  **
C*             NXA     NYA     SCRWL   SCRWR   SSX     T1      T2     **
C*             XSA     YSA     ZRATI   ZRATI2                         **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKIN  BLCKOT  FILSUF                                 **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines
      EXTERNAL BLCKIN , BLCKOT , FILSUF

C     Parameter defintions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF
      CHARACTER*78 COM(5)
      INTEGER LD , IRESRV , NLOAD , NLODF , KCOM , IPSYM , N1 , N2 , N , 
     &        NP , M1 , M2 , M , MP , IFAIL , DEBUG 
      INTEGER IP(2*LD) , ICON1(2*LD) , ICON2(2*LD) , 
     &        ITAG(2*LD) , ICONX(LD)
      REAL*8 EPSR , SIG , SCRWLT , SCRWRT , FMHZ , WLAM 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , BI(LD) , 
     &       ALP(LD) , BET(LD) , T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 CM(IRESRV) , ZARRAY(LD)

C     Local variables.
      INTEGER I , IOP , IOUT , J , K , NNEQ , NOP , NNPEQ

C     Common storage.
      INCLUDE 'gnd.inc'
      INCLUDE 'ggrid.inc'
      INCLUDE 'matpar.inc'
      INCLUDE 'smat.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING GFOUT'
      ENDIF

      IF ( DEBUG.GT.5) THEN
         WRITE(99,8999) LD, IRESRV, NLOAD, NLODF, KCOM, IPSYM, LD,
     &                  N1, N2, N, NP, M1, M2, M, MP , NRADL , KSYMP ,
     &                  IFAR , IPERF
         WRITE(99,FMT='(15(1X,I12))') ICASE , NBLOKS , NPBLK , NLAST ,
     &                                NBLSYM , NPSYM , NLSYM , IMAT , 
     &                                ICASX , NBBX , NPBX , NLBX , 
     &                                NBBL , NPBL , NLBL
         WRITE(99,9000) EPSR, SIG, SCRWLT, SCRWRT, FMHZ, WLAM , CL,CH ,
     &                  SCRWL , SCRWR , T2 , DBLE(ZRATI) , 
     &                  DIMAG(ZRATI) , DBLE(ZRATI2) , DIMAG(ZRATI2) ,
     &                  DBLE(FRATI) , DIMAG(FRATI) , DBLE(T1) , 
     &                  DIMAG(T1)
         DO 8001 I=1,LD
            WRITE(98,9001) ICONX(I) , SALP(I) , X(I) , Y(I) , Z(I) ,
     &                     SI(I) , BI(I) , ALP(I) , BET(I) , T2X(I) , 
     &                     T2Y(I) , T2Z(I) , DBLE(ZARRAY(I)) ,
     &                     DIMAG(ZARRAY(I))
 8001    CONTINUE
         DO 8002 I=1,2*LD
            WRITE(97,9002) IP(I), ICON1(I), ICON2(I) , ITAG(I)
 8002    CONTINUE
         DO 8003 I=1,IRESRV
            WRITE(96,9003) DBLE(CM(I)) , DIMAG(CM(I))
 8003    CONTINUE
         WRITE(95,9004) COM(1)
         WRITE(95,9004) COM(2)
         WRITE(95,9004) COM(3)
         WRITE(95,9004) COM(4)
         WRITE(95,9004) COM(5)
         DO 8005 I=1,16
            DO 8004 J=1,16
               WRITE(95,9005) DBLE(SSX(I,J)) , DIMAG(SSX(I,J))
 8004       CONTINUE
 8005    CONTINUE
 8999    FORMAT(19(1X,I12))
 9000    FORMAT(19(1X,E24.18))
 9001    FORMAT(1X,I12,13(1X,E24.18))
 9002    FORMAT(4(1X,I12))
 9003    FORMAT(2(1X,E24.18))
 9004    FORMAT(1X,'XXX',A78,'XXX')
 9005    FORMAT(2(1X,E24.18))
      ENDIF
 
      CALL FILSUF(JOBNAM,FILBUF,0,SUFNGF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      OPEN (UNIT=CHNGFL,FILE=FILBUF,FORM='UNFORMATTED',
     &      STATUS='NEW',ERR=1000)
      GOTO 1010
 1000 IFAIL = 48
      WRITE (CHRSLT,14) CHNGFL , FILBUF
      RETURN
 1010 CONTINUE

      NNEQ = N + 2*M
      NNPEQ = NP + 2*MP
      NOP = NNEQ/NNPEQ
      WRITE (CHNGFL) N , NP , M , MP , WLAM , FMHZ , IPSYM , KSYMP , 
     &             IPERF , NRADL , EPSR , SIG , SCRWLT , SCRWRT , 
     &             NLOAD , KCOM
      IF ( N.EQ.0 ) GOTO 1
      WRITE (CHNGFL) (X(I),I=1,N) , (Y(I),I=1,N) , (Z(I),I=1,N)
      WRITE (CHNGFL) (SI(I),I=1,N) , (BI(I),I=1,N) , (ALP(I),I=1,N)
      WRITE (CHNGFL) (BET(I),I=1,N) , (SALP(I),I=1,N)
      WRITE (CHNGFL) (ICON1(I),I=1,N) , (ICON2(I),I=1,N)
      WRITE (CHNGFL) (ITAG(I),I=1,N)
      IF ( NLOAD.GT.0 ) WRITE (CHNGFL) (ZARRAY(I),I=1,N)
    1 IF ( M.EQ.0 ) GOTO 2
      J = LD - M + 1
      WRITE (CHNGFL) (X(I),I=J,LD) , (Y(I),I=J,LD) , (Z(I),I=J,LD)
      WRITE (CHNGFL) (SI(I),I=J,LD) , (BI(I),I=J,LD) , (ALP(I),I=J,LD)
      WRITE (CHNGFL) (BET(I),I=J,LD) , (SALP(I),I=J,LD)

C     ERROR CORRECTED 11/20/89
 
      WRITE (CHNGFL) (T2X(I),I=J,LD) , (T2Y(I),I=J,LD)
      WRITE (CHNGFL) (T2Z(I),I=J,LD)
C      WRITE (CHNGFL) (ICON1(I),I=J,LD),(ICON2(I),I=J,LD)
C      WRITE (CHNGFL) (ITAG(I),I=J,LD)

    2 WRITE (CHNGFL) ICASE , NBLOKS , NPBLK , NLAST , NBLSYM , NPSYM , 
     &             NLSYM , IMAT
      IF ( IPERF.EQ.2 ) WRITE (CHNGFL) AR1 , AR2 , AR3 , EPSCF , DXA , 
     &                               DYA , XSA , YSA , NXA , NYA
      IF ( NOP.GT.1 ) WRITE (CHNGFL) ((SSX(I,J),I=1,NOP),J=1,NOP)
      WRITE (CHNGFL) (IP(I),I=1,NNEQ) , COM
      IF ( ICASE.GT.2 ) GOTO 3
      IOUT = NNEQ*NNPEQ
      WRITE (CHNGFL) (CM(I),I=1,IOUT)
      GOTO 12
    3 IF ( ICASE.NE.4 ) GOTO 5
      REWIND CHTMP3
      I = NNPEQ*NNPEQ
      DO 4 K = 1 , NOP
         READ (CHTMP3) (CM(J),J=1,I)
         WRITE (CHNGFL) (CM(J),J=1,I)
    4 CONTINUE
      REWIND CHTMP3
      GOTO 12
    5 REWIND CHTMP3
      REWIND CHTMP4
      IF ( ICASE.EQ.5 ) GOTO 8
      IOUT = NPBLK*NNEQ*2
      DO 6 I = 1 , NBLOKS
         CALL BLCKIN(CM,13,1,IOUT,1,201,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         CALL BLCKOT(CM,CHNGFL,1,IOUT,DEBUG)
    6 CONTINUE
      DO 7 I = 1 , NBLOKS
         CALL BLCKIN(CM,14,1,IOUT,1,203,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         CALL BLCKOT(CM,CHNGFL,1,IOUT,DEBUG)
    7 CONTINUE
      GOTO 12
    8 IOUT = NPSYM*NNPEQ*2
      DO 11 IOP = 1 , NOP
         DO 9 I = 1 , NBLSYM
            CALL BLCKIN(CM,13,1,IOUT,1,205,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            CALL BLCKOT(CM,CHNGFL,1,IOUT,DEBUG)
    9    CONTINUE
         DO 10 I = 1 , NBLSYM
            CALL BLCKIN(CM,14,1,IOUT,1,207,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            CALL BLCKOT(CM,CHNGFL,1,IOUT,DEBUG)
   10    CONTINUE
   11 CONTINUE
      REWIND CHTMP3
      REWIND CHTMP4
   12 REWIND CHNGFL
      WRITE (CHRSLT,13) CHNGFL , IMAT
 
      CLOSE (CHNGFL)
 
      RETURN

   13 FORMAT (///,' ****NUMERICAL GREEN''S FUNCTION FILE ON TAPE',I3,
     &        ' ****',/,5X,'MATRIX STORAGE -',I7,' COMPLEX NUMBERS',///)

   14 FORMAT (///,' ****NUMERICAL GREEN''S FUNCTION FILE ON TAPE',I3,
     &        ' **** FILE ',A,' ALREADY EXISTS',///)
 
      END
