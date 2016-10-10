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

      SUBROUTINE DATAGN( LD , SALP , IRESRV , CM , NLODF , ZARRAY , 
     &                   KCOM , COM , IP , EPSR , SIG , SCRWLT , 
     &                   SCRWRT , FMHZ , IPSYM , N1 , N2 ,N , NP , M1 , 
     &                   M2 , M , MP , ICON1 , ICON2 , ITAG , ICONX , 
     &                   WLAM , X , Y , Z , SI , BI , ALP , BET , X2 , 
     &                   Y2 , Z2 , CAB , SAB , T1X , T1Y , T1Z , T2X , 
     &                   T2Y , T2Z , IPLP1 , IPLP2 , JMAX , NSMAXX , 
     &                   NPMAX , NSCON , NPCON , JCO , ISCON , IPCON , 
     &                   JOBNAM , FILBUF , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* DATAGN is the main routine for input of geometry data.             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
C* INPUT  - IRESRV                                                    **
C* PASSED - CM                                                        **
C* PASSED - NLODF                                                     **
C* PASSED - ZARRAY                                                    **
C* PASSED - KCOM                                                      **
C* PASSED - COM                                                       **
C* PASSED - IP                                                        **
C* PASSED - EPSR                                                      **
C* PASSED - SIG                                                       **
C* PASSED - SCRWLT                                                    **
C* PASSED - SCRWRT                                                    **
C* PASSED - FMHZ                                                      **
C* OUTPUT - IPSYM                                                     **
C* OUTPUT - N1                                                        **
C* OUTPUT - N2                                                        **
C* OUTPUT - N                                                         **
C* OUTPUT - NP                                                        **
C* OUTPUT - M1                                                        **
C* OUTPUT - M2                                                        **
C* OUTPUT - M                                                         **
C* OUTPUT - MP                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - ITAG                                                      **
C* PASSED - ICONX                                                     **
C* PASSED - WLAM                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - SI                                                        **
C* OUTPUT - BI                                                        **
C* PASSED - ALP                                                       **
C* PASSED - BET                                                       **
C* OUTPUT - X2                                                        **
C* OUTPUT - Y2                                                        **
C* OUTPUT - Z2                                                        **
C* OUTPUT - CAB                                                       **
C* OUTPUT - SAB                                                       **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* OUTPUT - IPLP1                                                     **
C* OUTPUT - IPLP2                                                     **
C* INPUT  - JMAX                                                      **
C* INPUT  - NSMAXX                                                    **
C* INPUT  - NPMAX                                                     **
C* PASSED - NSCON                                                     **
C* PASSED - NPCON                                                     **
C* PASSED - JCO                                                       **
C* PASSED - ISCON                                                     **
C* PASSED - IPCON                                                     **
C* PASSED - JOBNAM                                                    **
C* PASSED - FILBUF                                                    **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ARC     ATGN2   CONECT  GFIL    HELIX   MOVE    PATCH  **
C*             READGM  REFLC   WIRE                                   **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      REAL*8 ATGN2
      EXTERNAL ARC , ATGN2 , CONECT , GFIL , HELIX , MOVE , PATCH ,
     &         READGM , REFLC , WIRE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF
      CHARACTER*78 COM(5)
      INTEGER LD , IRESRV , NLODF , KCOM , IPSYM , N1 , N2 , N , NP , 
     &        M1 , M2 , M , MP , IPLP1 , IPLP2 , JMAX , NSMAXX , NPMAX , 
     &        NSCON , NPCON , IFAIL , DEBUG
      INTEGER IP(2*LD) , ICON1(2*LD) , ICON2(2*LD) , 
     &        ITAG(2*LD) , ICONX(LD) , JCO(JMAX) , 
     &        ISCON(NSMAXX) , IPCON(NPMAX)
      REAL*8 EPSR , SIG , FMHZ , SCRWLT , SCRWRT , WLAM 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , ALP(LD) , BET(LD) , 
     &       X2(LD) , Y2(LD) , Z2(LD) , CAB(LD) , 
     &       SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , 
     &       T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 CM(IRESRV) , ZARRAY(LD)

C     Local variables.
      CHARACTER*2 GM
      CHARACTER*2 ATST(14)
      CHARACTER*1 IFX(2) , IFY(2) , IFZ(2) , IPT(4)
      INTEGER I , I1 , I2 , IPHD , IPSAV , ISCT , ITG , IX , IY , IZ , 
     &        J , MPSAV , NPSAV , NS , NWIRE
      REAL*8 DUMMY , RAD , TA , TD , X3 , X4 , XS1 , XXS2 , XW1 , XW2 , 
     &       Y3 , Y4 , YS1 , YS2 , YW1 , YW2 , Z3 , Z4 , ZS1 ,ZS2 , 
     &       ZW1 , ZW2
 
C     Data initialisation/
      DATA ATST/'GW' , 'GX' , 'GR' , 'GS' , 'GE' , 'GM' , 'SP' , 'SM' , 
     &          'GF' , 'GA' , 'SC' , 'GC' , 'GH' , 'IG'/
      DATA IFX/' ' , 'X'/
      DATA IFY/' ' , 'Y'/
      DATA IFZ/' ' , 'Z'/
      DATA TA/0.0174532925199432957692369D0/
      DATA TD/57.29577951308232087679815D0/
      DATA IPT/'P' , 'R' , 'T' , 'Q'/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING DATAGN'
      ENDIF
 
      IPSYM = 0
      NWIRE = 0
      N = 0
      NP = 0
      M = 0
      MP = 0
      N1 = 0
      N2 = 1
      M1 = 0
      M2 = 1
      ISCT = 0
      IPHD = 0
      X4 = 0.0D0
      Y4 = 0.0D0
      Z4 = 0.0D0
      X3 = 0.0D0
      Y3 = 0.0D0
      Z3 = 0.0D0

C     Read geometry data card and branch to section for operation
C     requested.

    1 CALL READGM(2,GM,ITG,NS,XW1,YW1,ZW1,XW2,YW2,ZW2,RAD,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
C IDF ignore comments.
      IF ( GM.EQ.ATST(14) ) GOTO 1
C end IDF.
      IF ( N+M.GT.LD ) GOTO 37
      IF ( GM.EQ.ATST(9) ) GOTO 27
      IF ( IPHD.EQ.1 ) GOTO 2
      WRITE (CHRSLT,40)
      WRITE (CHRSLT,41)
      IPHD = 1
    2 IF ( GM.EQ.ATST(11) ) GOTO 10
      ISCT = 0
      IF ( GM.EQ.ATST(1) ) GOTO 3
      IF ( GM.EQ.ATST(2) ) GOTO 18
      IF ( GM.EQ.ATST(3) ) GOTO 19
      IF ( GM.EQ.ATST(4) ) GOTO 21
      IF ( GM.EQ.ATST(7) ) GOTO 9
      IF ( GM.EQ.ATST(8) ) GOTO 13
      IF ( GM.EQ.ATST(5) ) GOTO 29
      IF ( GM.EQ.ATST(6) ) GOTO 26
      IF ( GM.EQ.ATST(10) ) GOTO 8
C***
      IF ( GM.EQ.ATST(13) ) GOTO 123
C***
      GOTO 36

C     Generate segment data for straight wire.

    3 NWIRE = NWIRE + 1
      I1 = N + 1
      I2 = N + NS
      WRITE (CHRSLT,43) NWIRE , XW1 , YW1 , ZW1 , XW2 , YW2 , ZW2 , 
     &                  RAD , NS , I1 , I2 , ITG
      IF ( RAD.EQ.0.0D0 ) GOTO 4
      XS1 = 1.0D0
      YS1 = 1.0D0
      GOTO 7
    4 CALL READGM(2,GM,IX,IY,XS1,YS1,ZS1,DUMMY,DUMMY,DUMMY,DUMMY,IFAIL,
     &            DEBUG)
      IF(IFAIL.NE.0) RETURN
C IDF ignore comments.
      IF ( GM.EQ.ATST(14) ) GOTO 4
C end IDF.

      IF ( GM.EQ.ATST(12) ) GOTO 6
    5 WRITE (CHRSLT,48)
      IFAIL=10
      RETURN
    6 WRITE (CHRSLT,61) XS1 , YS1 , ZS1
      IF ( YS1.EQ.0.0D0 .OR. ZS1.EQ.0.0D0 ) GOTO 5
      RAD = YS1
      YS1 = (ZS1/YS1)**(1.0D0/(NS-1.0D0))
    7 CALL WIRE(XW1,YW1,ZW1,XW2,YW2,ZW2,RAD,XS1,YS1,NS,ITG,LD,IPSYM,
     &          N,NP,M,MP,ITAG,X,Y,Z,BI,X2,Y2,Z2,DEBUG)
      GOTO 1

C     Generate segment data for wire arc.

    8 NWIRE = NWIRE + 1
      I1 = N + 1
      I2 = N + NS
      WRITE (CHRSLT,38) NWIRE , XW1 , YW1 , ZW1 , XW2 , NS , I1 , 
     &                  I2 , ITG
      CALL ARC(ITG,NS,XW1,YW1,ZW1,XW2,LD,IPSYM,N,NP,M,MP,ITAG,X,Y,
     &         Z,BI,X2,Y2,Z2,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 1

C     Generate helix.

  123 NWIRE = NWIRE + 1
      I1 = N + 1
      I2 = N + NS
      WRITE (CHRSLT,124) XW1 , YW1 , NWIRE , ZW1 , XW2 , YW2 , ZW2 , 
     &                   RAD , NS , I1 , I2 , ITG
      CALL HELIX(XW1,YW1,ZW1,XW2,YW2,ZW2,RAD,NS,ITG,LD,IPSYM,N,NP,M,
     &           MP,ITAG,X,Y,Z,BI,X2,Y2,Z2,DEBUG)
      GOTO 1

C     Generate single new patch.

    9 I1 = M + 1
      NS = NS + 1
      IF ( ITG.NE.0 ) GOTO 17
      WRITE (CHRSLT,51) I1 , IPT(NS) , XW1 , YW1 , ZW1 , XW2 , YW2 , ZW2
      IF ( NS.EQ.2 .OR. NS.EQ.4 ) ISCT = 1
      IF ( NS.GT.1 ) GOTO 14
      XW2 = XW2*TA
      YW2 = YW2*TA
      GOTO 16
   10 IF ( ISCT.EQ.0 ) GOTO 17
      I1 = M + 1
      NS = NS + 1
      IF ( ITG.NE.0 ) GOTO 17
      IF ( NS.NE.2 .AND. NS.NE.4 ) GOTO 17
      XS1 = X4
      YS1 = Y4
      ZS1 = Z4
      XXS2 = X3
      YS2 = Y3
      ZS2 = Z3
      X3 = XW1
      Y3 = YW1
      Z3 = ZW1
      IF ( NS.NE.4 ) GOTO 11
      X4 = XW2
      Y4 = YW2
      Z4 = ZW2
   11 XW1 = XS1
      YW1 = YS1
      ZW1 = ZS1
      XW2 = XXS2
      YW2 = YS2
      ZW2 = ZS2
      IF ( NS.EQ.4 ) GOTO 12
      X4 = XW1 + X3 - XW2
      Y4 = YW1 + Y3 - YW2
      Z4 = ZW1 + Z3 - ZW2
   12 WRITE (CHRSLT,51) I1 , IPT(NS) , XW1 , YW1 , ZW1 , XW2 , YW2 , ZW2
      WRITE (CHRSLT,39) X3 , Y3 , Z3 , X4 , Y4 , Z4
      GOTO 16

C     Generate multiple-patch surface.

   13 I1 = M + 1
      WRITE (CHRSLT,59) I1 , IPT(2) , XW1 , YW1 , ZW1 , XW2 , YW2 , 
     &                  ZW2 , ITG , NS
      IF ( ITG.LT.1 .OR. NS.LT.1 ) GOTO 17
   14 CALL READGM(2,GM,IX,IY,X3,Y3,Z3,X4,Y4,Z4,DUMMY,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
C IDF ignore comments.
      IF ( GM.EQ.ATST(14) ) GOTO 14
C end IDF.
      IF ( NS.NE.2 .AND. ITG.LT.1 ) GOTO 15
      X4 = XW1 + X3 - XW2
      Y4 = YW1 + Y3 - YW2
      Z4 = ZW1 + Z3 - ZW2
   15 WRITE (CHRSLT,39) X3 , Y3 , Z3 , X4 , Y4 , Z4
      IF ( GM.NE.ATST(11) ) GOTO 17
   16 CALL PATCH(ITG,NS,XW1,YW1,ZW1,XW2,YW2,ZW2,X3,Y3,Z3,X4,Y4,Z4,
     &           LD,SALP,IPSYM,N,NP,M,MP,X,Y,Z,BI,T1X,T1Y,T1Z,
     &           T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 1
   17 WRITE (CHRSLT,60)
      IFAIL=11
      RETURN

C     Reflect structure along X,Y, or Z axes or rotate to form cylinder.

   18 IY = NS/10
      IZ = NS - IY*10
      IX = IY/10
      IY = IY - IX*10
      IF ( IX.NE.0 ) IX = 1
      IF ( IY.NE.0 ) IY = 1
      IF ( IZ.NE.0 ) IZ = 1
      WRITE (CHRSLT,44) IFX(IX+1) , IFY(IY+1) , IFZ(IZ+1) , ITG
      GOTO 20
   19 WRITE (CHRSLT,45) NS , ITG
      IX = -1
   20 CALL REFLC(IX,IY,IZ,ITG,NS,LD,SALP,IPSYM,N1,N2,N,NP,M1,M2,M,MP,
     &           ITAG,X,Y,Z,BI,X2,Y2,Z2,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,
     &           DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 1

C     Scale structure dimensions by factor XW1.

   21 IF ( N.LT.N2 ) GOTO 23
      DO 22 I = N2 , N
         X(I) = X(I)*XW1
         Y(I) = Y(I)*XW1
         Z(I) = Z(I)*XW1
         X2(I) = X2(I)*XW1
         Y2(I) = Y2(I)*XW1
         Z2(I) = Z2(I)*XW1
         BI(I) = BI(I)*XW1
   22 CONTINUE
   23 IF ( M.LT.M2 ) GOTO 25
      YW1 = XW1*XW1
      IX = LD + 1 - M
      IY = LD - M1
      DO 24 I = IX , IY
         X(I) = X(I)*XW1
         Y(I) = Y(I)*XW1
         Z(I) = Z(I)*XW1
         BI(I) = BI(I)*YW1
   24 CONTINUE
   25 WRITE (CHRSLT,46) XW1
      GOTO 1

C     Move structure or reproduce original structure in new positions.

   26 WRITE (CHRSLT,47) ITG , NS , XW1 , YW1 , ZW1 , XW2 , YW2 , 
     &                  ZW2 , RAD
      XW1 = XW1*TA
      YW1 = YW1*TA
      ZW1 = ZW1*TA
      CALL MOVE(XW1,YW1,ZW1,XW2,YW2,ZW2,INT(RAD+0.5D0),NS,ITG,LD,
     &          SALP,IPSYM,N2,N,NP,M1,M2,M,MP,ITAG,X,Y,Z,BI,X2,Y2,Z2,
     &          T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 1

C     Read numerical green's function tape.

   27 IF ( N+M.EQ.0 ) GOTO 28
      WRITE (CHRSLT,52)
      IFAIL=12
      RETURN
   28 CALL GFIL(ITG,LD,SALP,IRESRV,CM,NLODF,ZARRAY,KCOM,COM,IP,EPSR,
     &          SIG,SCRWLT,SCRWRT,FMHZ,IPSYM,N1,N2,N,NP,M1,M2,M,MP,
     &          ICON1,ICON2,ITAG,WLAM,X,Y,Z,SI,BI,ALP,BET,T2X,T2Y,T2Z,
     &          JOBNAM,FILBUF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      NPSAV = NP
      MPSAV = MP
      IPSAV = IPSYM
      GOTO 1

C     Terminate structure geometry input.

   29 IF ( NS.EQ.0 ) GOTO 290
      IPLP1 = 1
      IPLP2 = 1
  290 IX = N1 + M1

      IF ( IX.EQ.0 ) GOTO 30
      NP = N
      MP = M
      IPSYM = 0
   30 CALL CONECT(ITG,LD,SALP,IPSYM,N1,N2,N,NP,M1,M2,M,MP,ICON1,ICON2,
     &            ICONX,X,Y,Z,BI,X2,Y2,Z2,T1X,T1Y,T1Z,T2X,T2Y,T2Z,JMAX,
     &            NSMAXX,NPMAX,NSCON,NPCON,JCO,ISCON,IPCON,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IF ( IX.EQ.0 ) GOTO 31
      NP = NPSAV
      MP = MPSAV
      IPSYM = IPSAV
   31 IF ( N+M.GT.LD ) GOTO 37
      IF ( N.EQ.0 ) GOTO 33
      WRITE (CHRSLT,53)
      WRITE (CHRSLT,54)
      DO 32 I = 1 , N
         XW1 = X2(I) - X(I)
         YW1 = Y2(I) - Y(I)
         ZW1 = Z2(I) - Z(I)
         X(I) = (X(I)+X2(I))*0.5D0
         Y(I) = (Y(I)+Y2(I))*0.5D0
         Z(I) = (Z(I)+Z2(I))*0.5D0
         XW2 = XW1*XW1 + YW1*YW1 + ZW1*ZW1
         YW2 = DSQRT(XW2)
         YW2 = (XW2/YW2+YW2)*0.5D0
         SI(I) = YW2
         CAB(I) = XW1/YW2
         SAB(I) = YW1/YW2
         XW2 = ZW1/YW2
         IF ( XW2.GT.1.0D0 ) XW2 = 1.0D0
         IF ( XW2.LT.-1.0D0 ) XW2 = -1.0D0
         SALP(I) = XW2
         XW2 = DASIN(XW2)*TD
         YW2 = ATGN2(YW1,XW1)*TD
         WRITE (CHRSLT,55) I , X(I) , Y(I) , Z(I) , SI(I) , XW2 , YW2 , 
     &                BI(I) , ICON1(I) , I , ICON2(I) , ITAG(I)

         IF ( IPLP1.NE.1 ) GOTO 320
         WRITE (CHRPAT,*) X(I) , Y(I) , Z(I) , SI(I) , XW2 , YW2 , 
     &                    BI(I) , ICON1(I) , I , ICON2(I)
  320    CONTINUE

         IF ( SI(I).GT.1.D-20 .AND. BI(I).GT.0. ) GOTO 32
         WRITE (CHRSLT,56)
         IFAIL=13
         RETURN
   32 CONTINUE
   33 IF ( M.EQ.0 ) GOTO 35
      WRITE (CHRSLT,57)
      J = LD + 1
      DO 34 I = 1 , M
         J = J - 1
         XW1 = (T1Y(J)*T2Z(J)-T1Z(J)*T2Y(J))*SALP(J)
         YW1 = (T1Z(J)*T2X(J)-T1X(J)*T2Z(J))*SALP(J)
         ZW1 = (T1X(J)*T2Y(J)-T1Y(J)*T2X(J))*SALP(J)
         WRITE (CHRSLT,58) I , X(J) , Y(J) , Z(J) , XW1 , YW1 , ZW1 , 
     &                     BI(J) , T1X(J) , T1Y(J) , T1Z(J) , T2X(J) , 
     &                     T2Y(J) , T2Z(J)
   34 CONTINUE
   35 RETURN
   36 WRITE (CHRSLT,48)
      WRITE (CHRSLT,49) GM , ITG , NS , XW1 , YW1 , ZW1 , XW2 , YW2 , 
     &                  ZW2 , RAD
      IFAIL=10
      RETURN
   37 WRITE (CHRSLT,50)
      IFAIL=14
      RETURN
      
  124 FORMAT (5X,'HELIX STRUCTURE-   AXIAL SPACING BETWEEN TURNS =',
     &        F8.3,' TOTAL AXIAL LENGTH =',F8.3/1X,I5,2X,
     &        'RADIUS OF HELIX =',4(2X,F8.3),7X,F11.5,I8,4X,I5,1X,I5,3X,
     &        I5)

   38 FORMAT (1X,I5,2X,'ARC RADIUS =',F9.5,2X,'FROM',F8.3,' TO',F8.3,
     &        ' DEGREES',11X,F11.5,2X,I5,4X,I5,1X,I5,3X,I5)
   39 FORMAT (6X,3F11.5,1X,3F11.5)
   40 FORMAT (////,33X,'- - - STRUCTURE SPECIFICATION - - -',//,37X,
     &        'COORDINATES MUST BE INPUT IN',/,37X,
     &        'METERS OR BE SCALED TO METERS',/,37X,
     &        'BEFORE STRUCTURE INPUT IS ENDED',//)
   41 FORMAT (2X,'WIRE',79X,'NO. OF',4X,'FIRST',2X,'LAST',5X,'TAG',/,2X,
     &        'NO.',8X,'X1',9X,'Y1',9X,'Z1',10X,'X2',9X,'Y2',9X,'Z2',6X,
     &        'RADIUS',3X,'SEG.',5X,'SEG.',3X,'SEG.',5X,'NO.')
C  42 FORMAT (A2,I3,I5,7F10.5)
   43 FORMAT (1X,I5,3F11.5,1X,4F11.5,2X,I5,4X,I5,1X,I5,3X,I5)
   44 FORMAT (6X,'STRUCTURE REFLECTED ALONG THE AXES',3(1X,A1),
     &        '.  TAGS INCREMENTED BY',I5)
   45 FORMAT (6X,'STRUCTURE ROTATED ABOUT Z-AXIS',I3,
     &        ' TIMES.  LABELS INCREMENTED BY',I5)
   46 FORMAT (6X,'STRUCTURE SCALED BY FACTOR',F10.5)
   47 FORMAT (6X,'THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS -'/6X,
     &        I3,I5,7F10.5)
   48 FORMAT (' GEOMETRY DATA CARD ERROR')
   49 FORMAT (1X,A2,I3,I5,7F10.5)
   50 FORMAT (
     &' NUMBER OF WIRE SEGMENTS AND SURFACE PATCHES EXCEEDS DIMENSION LI
     &MIT.')
   51 FORMAT (1X,I5,A1,F10.5,2F11.5,1X,3F11.5)
   52 FORMAT (' ERROR - GF MUST BE FIRST GEOMETRY DATA CARD')
   53 FORMAT (////33X,'- - - - SEGMENTATION DATA - - - -',//,40X,
     &        'COORDINATES IN METERS',//,25X,
     &        'I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I',//)
   54 FORMAT (2X,'SEG.',3X,'COORDINATES OF SEG. CENTER',5X,'SEG.',5X,
     &        'ORIENTATION ANGLES',4X,'WIRE',4X,'CONNECTION DATA',3X,
     &        'TAG',/,2X,'NO.',7X,'X',9X,'Y',9X,'Z',7X,'LENGTH',5X,
     &        'ALPHA',5X,'BETA',6X,'RADIUS',4X,'I-',3X,'I',4X,'I+',4X,
     &        'NO.')
   55 FORMAT (1X,I5,4F10.5,1X,3F10.5,1X,3I5,2X,I5)
   56 FORMAT (' SEGMENT DATA ERROR')
   57 FORMAT (////,44X,'- - - SURFACE PATCH DATA - - -',//,49X,
     &        'COORDINATES IN METERS',//,1X,'PATCH',5X,
     &        'COORD. OF PATCH CENTER',7X,'UNIT NORMAL VECTOR',6X,
     &        'PATCH',12X,'COMPONENTS OF UNIT TANGENT VECTORS',/,2X,
     &        'NO.',6X,'X',9X,'Y',9X,'Z',9X,'X',7X,'Y',7X,'Z',7X,'AREA',
     &        7X,'X1',6X,'Y1',6X,'Z1',7X,'X2',6X,'Y2',6X,'Z2')
   58 FORMAT (1X,I4,3F10.5,1X,3F8.4,F10.5,1X,3F8.4,1X,3F8.4)
   59 FORMAT (1X,I5,A1,F10.5,2F11.5,1X,3F11.5,5X,'SURFACE -',I4,' BY',
     &        I3,' PATCHES')
   60 FORMAT (' PATCH DATA ERROR')
   61 FORMAT (9X,'ABOVE WIRE IS TAPERED.  SEG. LENGTH RATIO =',F9.5,/,
     &        33X,'RADIUS FROM',F9.5,' TO',F9.5)
 
      END
