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

      SUBROUTINE RDPAT( LD , SALP , AIR , AII , BIR , BII , CIR , 
     &                  CII , CUR , SCRWLT , SCRWRT , N , M , WLAM , X , 
     &                  Y , Z , SI , GAIN , CAB , SAB , XS , YS , ZS , 
     &                  Z2S , IPLP1 , IPLP2 , IPLP3 , IPLP4 , NTH , 
     &                  NPH , IPD ,IAVP , INOR , IAX , IXTYP , THETS , 
     &                  PHIS ,DTH , DPH , RFLD , GNOR , CLT , CHT , 
     &                  EPSR2 , SIG2 , XPR6 , PINR , PNLR , PLOSS , 
     &                  JOBNAM , FILBUF , MULTIF , RDPCNT , IFAIL ,
     &                  DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute radiation pattern, gain, normalized gain.                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - LD                                                        **
C* PASSED - SALP                                                      **
C* PASSED - AIR                                                       **
C* PASSED - AII                                                       **
C* PASSED - BIR                                                       **
C* PASSED - BII                                                       **
C* PASSED - CIR                                                       **
C* PASSED - CII                                                       **
C* PASSED - CUR                                                       **
C* INPUT  - SCRWLT                                                    **
C* INPUT  - SCRWRT                                                    **
C* PASSED - N                                                         **
C* PASSED - M                                                         **
C* INPUT  - WLAM                                                      **
C* PASSED - X                                                         **
C* PASSED - Y                                                         **
C* PASSED - Z                                                         **
C* PASSED - SI                                                        **
C* OUTPUT - GAIN                                                      **
C* PASSED - CAB                                                       **
C* PASSED - SAB                                                       **
C* PASSED - XS                                                        **
C* PASSED - YS                                                        **
C* PASSED - ZS                                                        **
C* PASSED - Z2S                                                       **
C* INPUT  - IPLP1                                                     **
C* INPUT  - IPLP2                                                     **
C* INPUT  - IPLP3                                                     **
C* INPUT  - IPLP4                                                     **
C* INPUT  - NTH                                                       **
C* INPUT  - NPH                                                       **
C* INPUT  - IPD                                                       **
C* INPUT  - IAVP                                                      **
C* INPUT  - INOR                                                      **
C* INPUT  - IAX                                                       **
C* INPUT  - IXTYP                                                     **
C* INPUT  - THETS                                                     **
C* INPUT  - PHIS                                                      **
C* INPUT  - DTH                                                       **
C* INPUT  - DPH                                                       **
C* INPUT  - RFLD                                                      **
C* INPUT  - GNOR                                                      **
C* INPUT  - CLT                                                       **
C* INPUT  - CHT                                                       **
C* INPUT  - EPSR2                                                     **
C* INPUT  - SIG2                                                      **
C* INPUT  - XPR6                                                      **
C* OUTPUT - PINR                                                      **
C* INPUT  - PNLR                                                      **
C* INPUT  - PLOSS                                                     **
C* PASSED - JOBNAM                                                    **
C* INPUT  - FILBUF                                                    **
C* INPUT  - MULTIF                                                    **
C* OUTPUT - RDPCNT                                                    **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    CH      CL      ZRATI2                                 **
C* passes arg  KSYMP   ZRATI                                          **
C* uses value  IFAR    KSYMP   NRADL   ZRATI                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ATGN2   CANG    DB10    FFLD    FILSUF  GFLD           **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      REAL*8 ATGN2 , CANG , DB10
      EXTERNAL  ATGN2 , CANG , DB10 , FFLD , GFLD , FILSUF

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF
      INTEGER LD , N , M , IPLP1 , IPLP2 , IPLP3 , IPLP4 , 
     &        NTH , NPH , IPD , IAVP , INOR , IAX , IXTYP , MULTIF , 
     &        RDPCNT , IFAIL , DEBUG
      REAL*8 SCRWLT , SCRWRT , WLAM , THETS , PHIS ,DTH , DPH , RFLD , 
     &       GNOR , CLT , CHT , EPSR2 , SIG2 , XPR6 , PINR , PNLR ,
     &       PLOSS
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       GAIN(4*LD) , CAB(LD) , SAB(LD) , XS(LD) , YS(LD) , 
     &       ZS(LD) , Z2S(LD) 
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      CHARACTER*6 IGNTP(10) , IGAX(4) , IGTP(4) ,HPOL(3) , HBLK , HCIR , 
     &            HCLIF , ISENS
      INTEGER I , ITMP1 , ITMP2 , ITMP3 , ITMP4 , J , KPH , KTH , 
     &        NORMAX
      REAL*8 AXRAT , CDFAZ , DA , DFAZ , DFAZ2 , EMAJR2 , EMINR2 , 
     &       EPHA , EPHM , EPHM2 , ERDA , ERDM , ETHA , ETHM , ETHM2 , 
     &       EXRA , EXRM , GCON , GCOP , GMAX , GNH , GNMJ , GNMN , 
     &       GNV , GTOT , PHA , PHI , PI , PINT , PRAD , STILTA , TA , 
     &       TD , THA , THET , TILTA , TMP1 , TMP2 , TMP3 , TMP4 , 
     &       TMP5 , TMP6 , TSTOR1 , TSTOR2
      COMPLEX*16 ETH , EPH , ERD

C     Common storage.
      INCLUDE 'gnd.inc'
      
C     Data initialisation.
      DATA HPOL /'LINEAR' , 'RIGHT' , 'LEFT'/
      DATA HBLK /' '/
      DATA HCIR /'CIRCLE'/
      DATA IGTP /'    - ' , 'POWER ' , '- DIRE' , 'CTIVE '/
      DATA IGAX /' MAJOR' , ' MINOR' , ' VERT.' , ' HOR. '/
      DATA IGNTP /' MAJOR' , ' AXIS ' , ' MINOR' , ' AXIS ' , '   VER' , 
     &            'TICAL ' , ' HORIZ' , 'ONTAL ' , '      ' , 'TOTAL '/
      DATA PI /3.141592653589793238462643D0/
      DATA TA /0.0174532925199432957692369D0/
      DATA TD /57.29577951308232087679815D0/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING RDPAT'
      ENDIF

      IF ( MULTIF.EQ.1 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,RDPCNT,SUFRDP,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         OPEN(UNIT=CHRDPT,STATUS='NEW',FILE=FILBUF)
         RDPCNT = RDPCNT + 1
      ELSE IF ( MULTIF.EQ.2 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,0,SUFRDP,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF(RDPCNT.EQ.1) THEN
            OPEN(UNIT=CHRDPT,STATUS='NEW',FILE=FILBUF)
         ENDIF
         RDPCNT = RDPCNT + 1
      ENDIF
      NORMAX = 4 * LD
      IF ( IFAR.LT.2 ) GOTO 2
      WRITE (CHRSLT,36)
      IF ( IFAR.LE.3 ) GOTO 1
      WRITE (CHRSLT,37) NRADL , SCRWLT , SCRWRT
      IF ( IFAR.EQ.4 ) GOTO 2
    1 IF ( IFAR.EQ.2 .OR. IFAR.EQ.5 ) HCLIF = HPOL(1)
      IF ( IFAR.EQ.3 .OR. IFAR.EQ.6 ) HCLIF = HCIR
      CL = CLT/WLAM
      CH = CHT/WLAM
      ZRATI2 = CDSQRT(1.0D0/DCMPLX(EPSR2,-SIG2*WLAM*59.9584916D0))
      WRITE (CHRSLT,38) HCLIF , CLT , CHT , EPSR2 , SIG2
    2 IF ( IFAR.NE.1 ) GOTO 3
      WRITE (CHRSLT,42)
      GOTO 5
    3 I = 2*IPD + 1
      J = I + 1
      ITMP1 = 2*IAX + 1
      ITMP2 = ITMP1 + 1
      WRITE (CHRSLT,39)
      IF ( RFLD.LT.1.D-20 ) GOTO 4
      EXRM = 1.0D0/RFLD
      EXRA = RFLD/WLAM
      EXRA = -360.0D0*(EXRA-AINT(EXRA))
      WRITE (CHRSLT,40) RFLD , EXRM , EXRA
    4 WRITE (CHRSLT,41) IGTP(I) , IGTP(J) , IGAX(ITMP1) , IGAX(ITMP2)
    5 IF ( IXTYP.EQ.0 .OR. IXTYP.EQ.5 ) GOTO 7
      IF ( IXTYP.EQ.4 ) GOTO 6
      PRAD = 0.0D0
      GCON = 4.0D0*PI/(1.0D0+XPR6*XPR6)
      GCOP = GCON
      GOTO 8
    6 PINR = 394.5110617186928930645169D0*XPR6*XPR6*WLAM*WLAM
    7 GCOP = WLAM*WLAM*2.0D0*PI/(376.7303134617706554681984D0*PINR)
      PRAD = PINR - PLOSS - PNLR
      GCON = GCOP
      IF ( IPD.NE.0 ) GCON = GCON*PINR/PRAD
    8 I = 0
      GMAX = -1.0D10
      PINT = 0.0D0
      TMP1 = DPH*TA
      TMP2 = 0.5D0*DTH*TA
      PHI = PHIS - DPH
      DO 30 KPH = 1 , NPH
         PHI = PHI + DPH
         PHA = PHI*TA
         THET = THETS - DTH
         DO 29 KTH = 1 , NTH
            THET = THET + DTH
            IF ( KSYMP.EQ.2 .AND. THET.GT.90.00001D0 .AND. IFAR.NE.1 )
     &           GOTO 29
            THA = THET*TA
            IF ( IFAR.EQ.1 ) GOTO 9
            CALL FFLD(THA,PHA,ETH,EPH,LD,SALP,AIR,AII,BIR,BII,CIR, 
     &                CII,CUR,N,M,X,Y,Z,SI,CAB,SAB,XS,YS,ZS,Z2S)
            GOTO 10
    9       CALL GFLD(RFLD/WLAM,PHA,THET/WLAM,ETH,EPH,ERD,ZRATI,KSYMP,
     &                LD,SALP,AIR,AII,BIR,BII,CIR,CII,CUR,N,M,X,
     &                Y,Z,SI,CAB,SAB,XS,YS,ZS,Z2S)
            ERDM = CDABS(ERD)
            ERDA = CANG(ERD)
   10       ETHM2 = DBLE(ETH*DCONJG(ETH))
            ETHM = DSQRT(ETHM2)
            ETHA = CANG(ETH)
            EPHM2 = DBLE(EPH*DCONJG(EPH))
            EPHM = DSQRT(EPHM2)
            EPHA = CANG(EPH)
            IF ( IFAR.EQ.1 ) GOTO 28
C     Elliptical polarization calc.
            IF ( ETHM2.GT.1.D-20 .OR. EPHM2.GT.1.D-20 ) GOTO 11
            TILTA = 0.0D0
            EMAJR2 = 0.0D0
            EMINR2 = 0.0D0
            AXRAT = 0.0D0
            ISENS = HBLK
            GOTO 16
   11       DFAZ = EPHA - ETHA
            IF ( EPHA.LT.0. ) GOTO 12
            DFAZ2 = DFAZ - 360.0D0
            GOTO 13
   12       DFAZ2 = DFAZ + 360.0D0
   13       IF ( DABS(DFAZ).GT.DABS(DFAZ2) ) DFAZ = DFAZ2
            CDFAZ = DCOS(DFAZ*TA)
            TSTOR1 = ETHM2 - EPHM2
            TSTOR2 = 2.0D0*EPHM*ETHM*CDFAZ
            TILTA = 0.5D0*ATGN2(TSTOR2,TSTOR1)
            STILTA = DSIN(TILTA)
            TSTOR1 = TSTOR1*STILTA*STILTA
            TSTOR2 = TSTOR2*STILTA*DCOS(TILTA)
            EMAJR2 = -TSTOR1 + TSTOR2 + ETHM2
            EMINR2 = TSTOR1 - TSTOR2 + EPHM2
            IF ( EMINR2.LT.0. ) EMINR2 = 0.0D0
            AXRAT = DSQRT(EMINR2/EMAJR2)
            TILTA = TILTA*TD
            IF ( AXRAT.GT.1.D-5 ) GOTO 14
            ISENS = HPOL(1)
            GOTO 16
   14       IF ( DFAZ.GT.0. ) GOTO 15
            ISENS = HPOL(2)
            GOTO 16
   15       ISENS = HPOL(3)
   16       GNMJ = DB10(GCON*EMAJR2)
            GNMN = DB10(GCON*EMINR2)
            GNV = DB10(GCON*ETHM2)
            GNH = DB10(GCON*EPHM2)
            GTOT = DB10(GCON*(ETHM2+EPHM2))
            IF ( INOR.LT.1 ) GOTO 23
            I = I + 1
            IF ( I.GT.NORMAX ) GOTO 23
            GOTO (17,18,19,20,21) , INOR
   17       TSTOR1 = GNMJ
            GOTO 22
   18       TSTOR1 = GNMN
            GOTO 22
   19       TSTOR1 = GNV
            GOTO 22
   20       TSTOR1 = GNH
            GOTO 22
   21       TSTOR1 = GTOT
   22       GAIN(I) = TSTOR1
            IF ( TSTOR1.GT.GMAX ) GMAX = TSTOR1
   23       IF ( IAVP.EQ.0 ) GOTO 24
            TSTOR1 = GCOP*(ETHM2+EPHM2)
            TMP3 = THA - TMP2
            TMP4 = THA + TMP2
            IF ( KTH.EQ.1 ) TMP3 = THA
            IF ( KTH.EQ.NTH ) TMP4 = THA
            DA = DABS(TMP1*(DCOS(TMP3)-DCOS(TMP4)))
            IF ( KPH.EQ.1 .OR. KPH.EQ.NPH ) DA = 0.5D0*DA
            PINT = PINT + TSTOR1*DA
            IF ( IAVP.EQ.2 ) GOTO 29
   24       IF ( IAX.EQ.1 ) GOTO 25
            TMP5 = GNMJ
            TMP6 = GNMN
            GOTO 26
   25       TMP5 = GNV
            TMP6 = GNH
   26       ETHM = ETHM*WLAM
            EPHM = EPHM*WLAM
            IF ( RFLD.LT.1.D-20 ) GOTO 27
            ETHM = ETHM*EXRM
            ETHA = ETHA + EXRA
            EPHM = EPHM*EXRM
            EPHA = EPHA + EXRA
   27       WRITE (CHRSLT,43) THET , PHI , TMP5 , TMP6 , GTOT , AXRAT , 
     &                        TILTA , ISENS , ETHM , ETHA , EPHM , EPHA
            IF ( MULTIF.GT.0 ) THEN
               WRITE (CHRDPT,48) THET , PHI , TMP5 , TMP6 , GTOT , 
     &                           AXRAT , TILTA , ETHM , ETHA , 
     &                           EPHM , EPHA
            ENDIF
            
C      GO TO 29
C28    WRITE(CHRSLT,43)  RFLD,PHI,THET,ETHM,ETHA,EPHM,EPHA,ERDM,ERDA
            IF ( IPLP1.NE.3 ) GOTO 299
            IF ( IPLP3.EQ.0 ) GOTO 290
            IF ( IPLP2.EQ.1 .AND. IPLP3.EQ.1 ) WRITE (CHRPAT,*) THET , 
     &                                                      ETHM , ETHA
            IF ( IPLP2.EQ.1 .AND. IPLP3.EQ.2 ) WRITE (CHRPAT,*) THET , 
     &                                                      EPHM , EPHA
            IF ( IPLP2.EQ.2 .AND. IPLP3.EQ.1 ) WRITE (CHRPAT,*) PHI , 
     &                                                      ETHM , ETHA
            IF ( IPLP2.EQ.2 .AND. IPLP3.EQ.2 ) WRITE (CHRPAT,*) PHI , 
     &                                                      EPHM , EPHA
            IF ( IPLP4.EQ.0 ) GOTO 299
  290       IF ( IPLP2.EQ.1 .AND. IPLP4.EQ.1 ) WRITE (CHRPAT,*) THET , 
     &                                                          TMP5
            IF ( IPLP2.EQ.1 .AND. IPLP4.EQ.2 ) WRITE (CHRPAT,*) THET , 
     &                                                          TMP6
            IF ( IPLP2.EQ.1 .AND. IPLP4.EQ.3 ) WRITE (CHRPAT,*) THET , 
     &                                                          GTOT
            IF ( IPLP2.EQ.1 .AND. IPLP4.EQ.4 ) WRITE (CHRPAT,*) THET , 
     &                                            TMP5 , TMP6 , GTOT
            IF ( IPLP2.EQ.2 .AND. IPLP4.EQ.1 ) WRITE (CHRPAT,*) PHI , 
     &                                                          TMP5
            IF ( IPLP2.EQ.2 .AND. IPLP4.EQ.2 ) WRITE (CHRPAT,*) PHI , 
     &                                                          TMP6
            IF ( IPLP2.EQ.2 .AND. IPLP4.EQ.3 ) WRITE (CHRPAT,*) PHI , 
     &                                                          GTOT
            IF ( IPLP2.EQ.2 .AND. IPLP4.EQ.4 ) WRITE (CHRPAT,*) PHI , 
     &                                            TMP5 , TMP6 , GTOT
            GOTO 299
 28         WRITE (CHRSLT,44) RFLD , PHI , THET , ETHM , ETHA , EPHM , 
     &                        EPHA , ERDM , ERDA
            IF ( MULTIF.GT.0 ) THEN
               WRITE (CHRDPT,44) RFLD , PHI , THET , ETHM , ETHA , 
     &                           EPHM , EPHA , ERDM , ERDA
            ENDIF
  299       CONTINUE

   29    CONTINUE
   30 CONTINUE
      IF ( IAVP.EQ.0 ) GOTO 31
      TMP3 = THETS*TA
      TMP4 = TMP3 + DTH*TA*DBLE(NTH-1)
      TMP3 = DABS(DPH*TA*DBLE(NPH-1)*(DCOS(TMP3)-DCOS(TMP4)))
      PINT = PINT/TMP3
      TMP3 = TMP3/PI
      WRITE (CHRSLT,45) PINT , TMP3
   31 IF ( INOR.EQ.0 ) GOTO 35
      IF ( DABS(GNOR).GT.1.D-20 ) GMAX = GNOR
      ITMP1 = (INOR-1)*2 + 1
      ITMP2 = ITMP1 + 1
      WRITE (CHRSLT,46) IGNTP(ITMP1) , IGNTP(ITMP2) , GMAX
      ITMP2 = NPH*NTH
      IF ( ITMP2.GT.NORMAX ) ITMP2 = NORMAX
      ITMP1 = (ITMP2+2)/3
      ITMP2 = ITMP1*3 - ITMP2
      ITMP3 = ITMP1
      ITMP4 = 2*ITMP1
      IF ( ITMP2.EQ.2 ) ITMP4 = ITMP4 - 1
      DO 32 I = 1 , ITMP1
         ITMP3 = ITMP3 + 1
         ITMP4 = ITMP4 + 1
         J = (I-1)/NTH
         TMP1 = THETS + DBLE(I-J*NTH-1)*DTH
         TMP2 = PHIS + DBLE(J)*DPH
         J = (ITMP3-1)/NTH
         TMP3 = THETS + DBLE(ITMP3-J*NTH-1)*DTH
         TMP4 = PHIS + DBLE(J)*DPH
         J = (ITMP4-1)/NTH
         TMP5 = THETS + DBLE(ITMP4-J*NTH-1)*DTH
         TMP6 = PHIS + DBLE(J)*DPH
         TSTOR1 = GAIN(I) - GMAX
         IF ( I.EQ.ITMP1 .AND. ITMP2.NE.0 ) GOTO 33
         TSTOR2 = GAIN(ITMP3) - GMAX
         PINT = GAIN(ITMP4) - GMAX
         WRITE (CHRSLT,47) TMP1 , TMP2 , TSTOR1 , TMP3 , TMP4 , TSTOR2 , 
     &                TMP5 , TMP6 , PINT
   32 CONTINUE
      GOTO 35
   33 IF ( ITMP2.EQ.2 ) GOTO 34
      TSTOR2 = GAIN(ITMP3) - GMAX
      WRITE (CHRSLT,47) TMP1 , TMP2 , TSTOR1 , TMP3 , TMP4 , TSTOR2
      GOTO 35
   34 WRITE (CHRSLT,47) TMP1 , TMP2 , TSTOR1

      IF ( MULTIF.EQ.1 ) THEN
         CLOSE(CHRDPT)
      ENDIF
      
   35 RETURN
 
   36 FORMAT (///,31X,'- - - FAR FIELD GROUND PARAMETERS - - -',//)
   37 FORMAT (40X,'RADIAL WIRE GROUND SCREEN',/,40X,I5,' WIRES',/,40X,
     &        'WIRE LENGTH=',F8.2,' METERS',/,40X,'WIRE RADIUS=',1P,
     &        E10.3,' METERS')
   38 FORMAT (40X,A6,' CLIFF',/,40X,'EDGE DISTANCE=',F9.2,' METERS',/,
     &        40X,'HEIGHT=',F8.2,' METERS',/,40X,'SECOND MEDIUM -',/,
     &        40X,'RELATIVE DIELECTRIC CONST.=',F7.3,/,40X,
     &        'CONDUCTIVITY=',1P,E10.3,' MHOS')
   39 FORMAT (///,48X,'- - - RADIATION PATTERNS - - -')
   40 FORMAT (54X,'RANGE=',1P,E13.6,' METERS',/,54X,'EXP(-JKR)/R=',
     &        E12.5,' AT PHASE',0P,F7.2,' DEGREES',/)
   41 FORMAT (/,2X,'- - ANGLES - -',7X,2A6,'GAINS -',7X,
     &        '- - - POLARIZATION - - -',4X,'- - - E(THETA) - - -',4X,
     &        '- - - E(PHI) - - -',/,2X,'THETA',5X,'PHI',7X,A6,2X,A6,3X,
     &        'TOTAL',6X,'AXIAL',5X,'TILT',3X,'SENSE',
     &        2(5X,'MAGNITUDE',4X,'PHASE '),/,2(1X,'DEGREES',1X),
     &        3(6X,'DB'),8X,'RATIO',5X,'DEG.',8X,
     &        2(6X,'VOLTS/M',4X,'DEGREES'))
   42 FORMAT (///,28X,' - - - RADIATED FIELDS NEAR GROUND - - -',//,8X,
     &        '- - - LOCATION - - -',10X,'- - E(THETA) - -',8X,
     &        '- - E(PHI) - -',8X,'- - E(RADIAL) - -',/,7X,'RHO',6X,
     &        'PHI',9X,'Z',12X,'MAG',6X,'PHASE',9X,'MAG',6X,'PHASE',9X,
     &        'MAG',6X,'PHASE',/,5X,'METERS',3X,'DEGREES',4X,'METERS',
     &        8X,'VOLTS/M',3X,'DEGREES',6X,'VOLTS/M',3X,'DEGREES',6X,
     &        'VOLTS/M',3X,'DEGREES',/)
   43 FORMAT (1X,F7.2,F9.2,3X,3F8.2,F11.5,F9.2,2X,A6,2(1P,E15.5,0P,F9.2)
     &        )
   44 FORMAT (3X,F9.2,2X,F7.2,2X,F9.2,1X,3(3X,1P,E11.4,2X,0P,F7.2))
   45 FORMAT (//,3X,'AVERAGE POWER GAIN=',1P,E12.5,7X,
     &        'SOLID ANGLE USED IN AVERAGING=(',0P,F7.4,
     &        ')*PI STERADIANS.',//)
   46 FORMAT (//,37X,'- - - - NORMALIZED GAIN - - - -',//,37X,2A6,
     &        'GAIN',/,38X,'NORMALIZATION FACTOR =',F9.2,' DB',//,
     &        3(4X,'- - ANGLES - -',6X,'GAIN',7X),/,
     &        3(4X,'THETA',5X,'PHI',8X,'DB',8X),/,
     &        3(3X,'DEGREES',2X,'DEGREES',16X))
   47 FORMAT (3(1X,2F9.2,1X,F9.2,6X))
   48 FORMAT (1X,F7.2,F9.2,3X,3F8.2,F11.5,F9.2,2X,2(1P,E15.5,0P,F9.2))
      END
