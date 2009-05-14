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
C*                                                                    **
C* nec2d.F - Double precision NEC2 code as modified by I D Flintoft   **
C*                                                                    **
C* For authors and copyright information see the file COPYING in the  **
C* top directory of this distribution.                                **
C*                                                                    **

#include <config.h>

C*--------------------------------------------------------------------**

      SUBROUTINE NEC2D( JOBNAM, FILBUF, IFAIL , DEBUG , MULTIF , LD , 
     &                  IRESRV , NSMAX , LOADMX , NETMX , NORMF , 
     &                  MXCOUP , MXCOU2 , JMAX , NSMAXX ,  NPMAX , IP , 
     &                  ICON1 , ICON2 , ITAG , ICONX , IVQD , ISANT , 
     &                  IQDS , NCTAG , NCSEG , ISEG1 , ISEG2 , NTYP , 
     &                  JCO , ISCON , IPCON ,LDTYP , LDTAG , LDTAGF , 
     &                  LDTAGT , IX , SALP , AIR , AII , BIR , BII , 
     &                  CIR , CII , GAIN , X , Y  , Z , SI , BI , ALP , 
     &                  BET , X2 , Y2 , Z2 , CAB , SAB , T1X , T1Y , 
     &                  T1Z , T2X , T2Y , T2Z , XS , YS , ZS , Z2S , 
     &                  X11R , X11I , X12R , X12I , X22R , X22I , AX , 
     &                  BX , CX , ZLR , ZLI , ZLC , FNORM , XTEMP , 
     &                  YTEMP , ZTEMP , SITEMP , BITEMP , CM , CUR , 
     &                  ZARRAY , D , VQD , VSANT , VQDS , Y11A , Y12A , 
     &                  IPNT , NTEQA , NTSCA , CMN , RHNT , RHS , VSRC , 
     &                  RHNX )

C*--------------------------------------------------------------------**
C*                                                                    **
C* NUMERICAL ELECTROMAGNETICS CODE (NEC2)  DEVELOPED AT LAWRENCE      **
C* LIVERMORE LAB., LIVERMORE, CA.  (CONTACT G. BURKE AT 415-422-8414  **
C* FOR PROBLEMS WITH THE NEC CODE.)                                   **
C* FILE CREATED 4/11/80.                                              **
C*                                                                    **
C*                   ***********NOTICE**********                      **
C* THIS COMPUTER CODE MATERIAL WAS PREPARED AS AN ACCOUNT OF WORK     **
C* SPONSORED BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED     **
C* STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF      **
C* THEIR EMPLOYEES, NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR  **
C* THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR        **
C* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,    **
C* COMPLETENESS OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT  **
C* OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT         **
C* INFRINGE PRIVATELY-OWNED RIGHTS.                                   **
C*                                                                    **
C* DOUBLE PRECISION 6/4/85                                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - JOBNAM                                                    **
C* INPUT  - FILBUF                                                    **
C* OUTPUT - IFAIL                                                     **
C* INPUT  - MULTIF                                                    **
C* INPUT  - LD                                                        **
C* INPUT  - IRESRV                                                    **
C* INPUT  - NSMAX                                                     **
C* INPUT  - LOADMX                                                    **
C* INPUT  - NETMX                                                     **
C* INPUT  - NORMF                                                     **
C* INPUT  - MXCOUP                                                    **
C* INPUT  - MXCOU2                                                    **
C* INPUT  - JMAX                                                      **
C* INPUT  - NSMAXX                                                    **
C* INPUT  - NPMAX                                                     **
C* PASSED - IP                                                        **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - ITAG                                                      **
C* PASSED - ICONX                                                     **
C* OUTPUT - IVQD                                                      **
C* OUTPUT - ISANT                                                     **
C* PASSED - IQDS                                                      **
C* OUTPUT - NCTAG                                                     **
C* OUTPUT - NCSEG                                                     **
C* OUTPUT - ISEG1                                                     **
C* OUTPUT - ISEG2                                                     **
C* OUTPUT - NTYP                                                      **
C* PASSED - JCO                                                       **
C* PASSED - ISCON                                                     **
C* PASSED - IPCON                                                     **
C* OUTPUT - LDTYP                                                     **
C* OUTPUT - LDTAG                                                     **
C* OUTPUT - LDTAGF                                                    **
C* OUTPUT - LDTAGT                                                    **
C* PASSED - IX                                                        **
C* PASSED - SALP                                                      **
C* PASSED - AIR                                                       **
C* PASSED - AII                                                       **
C* INPUT  - BIR                                                       **
C* INPUT  - BII                                                       **
C* PASSED - CIR                                                       **
C* PASSED - CII                                                       **
C* PASSED - GAIN                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - SI                                                        **
C* OUTPUT - BI                                                        **
C* PASSED - ALP                                                       **
C* PASSED - BET                                                       **
C* PASSED - X2                                                        **
C* PASSED - Y2                                                        **
C* PASSED - Z2                                                        **
C* PASSED - CAB                                                       **
C* PASSED - SAB                                                       **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* PASSED - XS                                                        **
C* PASSED - YS                                                        **
C* PASSED - ZS                                                        **
C* PASSED - Z2S                                                       **
C* OUTPUT - X11R                                                      **
C* OUTPUT - X11I                                                      **
C* OUTPUT - X12R                                                      **
C* OUTPUT - X12I                                                      **
C* OUTPUT - X22R                                                      **
C* OUTPUT - X22I                                                      **
C* PASSED - AX                                                        **
C* PASSED - BX                                                        **
C* PASSED - CX                                                        **
C* OUTPUT - ZLR                                                       **
C* OUTPUT - ZLI                                                       **
C* OUTPUT - ZLC                                                       **
C* OUTPUT - FNORM                                                     **
C* OUTPUT - XTEMP                                                     **
C* OUTPUT - YTEMP                                                     **
C* OUTPUT - ZTEMP                                                     **
C* OUTPUT - SITEMP                                                    **
C* OUTPUT - BITEMP                                                    **
C* PASSED - CM                                                        **
C* INPUT  - CUR                                                       **
C* OUTPUT - ZARRAY                                                    **
C* PASSED - D                                                         **
C* OUTPUT - VQD                                                       **
C* OUTPUT - VSANT                                                     **
C* PASSED - VQDS                                                      **
C* PASSED - Y11A                                                      **
C* PASSED - Y12A                                                      **
C* PASSED - IPNT                                                      **
C* PASSED - NTEQA                                                     **
C* PASSED - NTSCA                                                     **
C* PASSED - CMN                                                       **
C* PASSED - RHNT                                                      **
C* PASSED - RHS                                                       **
C* PASSED - VSRC                                                      **
C* PASSED - RHNX                                                      **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    AR1     AR2     AR3     DXA     DYA     EPSCF   FRATI  **
C*             ICASX   IFAR    IMAT    IPERF   KSYMP   NRADL   NXA    **
C*             NYA     SCRWL   SCRWR   T1      T2      U       U2     **
C*             XSA     YSA     ZRATI                                  **
C* passes arg  NPBX                                                   **
C* uses value  EPSCF   ICASX   IFAR    IMAT    IPERF   KSYMP   NPBX   **
C*             NRADL   NXA     SCRWR   U       ZRATI                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       CANG    CMNGF   CMSET   COUPLE  DATAGN  DB20    ERROR  **
C*             ETMNS   FACGF   FACTRS  FBLOCK  FBNGF   FILSUF  GFOUT  **
C*             ISEGNO  LOAD    NETWK   NFPAT   RDPAT   READMN  TIMER  **
C*             UPCASE                                                 **
C* called by   FDRIVE                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      REAL*8 CANG , DB20
      EXTERNAL CANG , CMNGF , CMSET , COUPLE , DATAGN , DB20 , ETMNS , 
     &         FACGF , FACTRS , FBLOCK , FBNGF , GFOUT , ISEGNO , LOAD , 
     &         NETWK , NFPAT , RDPAT , READMN , TIMER , UPCASE , FILSUF

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM, FILBUF
      INTEGER LD , IRESRV , NSMAX , LOADMX , NETMX , NORMF ,  
     &        MXCOUP , MXCOU2 , JMAX , NSMAXX ,  NPMAX , MULTIF , 
     &        IFAIL , NRECNT , NRHCNT , RDPCNT , DEBUG
      INTEGER IP(2*LD) , ICON1(2*LD) , ICON2(2*LD) , ITAG(2*LD) , 
     &        ICONX(LD) , IVQD(NSMAX) , ISANT(NSMAX) , IQDS(NSMAX) , 
     &        NCTAG(MXCOUP) , NCSEG(MXCOUP) , ISEG1(NETMX) , 
     &        ISEG2(NETMX) , NTYP(NETMX) , JCO(JMAX) , ISCON(NSMAXX) , 
     &        IPCON(NPMAX) , LDTYP(LOADMX) , LDTAG(LOADMX) , 
     &        LDTAGF(LOADMX) , LDTAGT(LOADMX) , IX(2*LD) , IPNT(NETMX) , 
     &        NTEQA(NETMX) , NTSCA(NETMX)
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , GAIN(4*LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , ALP(LD) , BET(LD) , X2(LD) , Y2(LD) , 
     &       Z2(LD) , CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , 
     &       T2X(LD) , T2Y(LD) , T2Z(LD) , XS(LD) , YS(LD) , ZS(LD) , 
     &       Z2S(LD) , X11R(NETMX) , X11I(NETMX) , X12R(NETMX) , 
     &       X12I(NETMX) , X22R(NETMX) , X22I(NETMX) , AX(JMAX) , 
     &       BX(JMAX) , CX(JMAX) , ZLR(LOADMX) , ZLI(LOADMX) , 
     &       ZLC(LOADMX) , FNORM(NORMF) , XTEMP(LD) , YTEMP(LD) , 
     &       ZTEMP(LD) , SITEMP(LD) , BITEMP(LD)
      COMPLEX*16 CM(IRESRV) , CUR(3*LD) , ZARRAY(LD) , D(2*LD) , 
     &           VQD(NSMAX) , VSANT(NSMAX) , VQDS(NSMAX) , 
     &           Y11A(MXCOU2) , Y12A(MXCOU2) , CMN(NETMX,NETMX) , 
     &           RHNT(NETMX) , RHS(3*LD) , VSRC(NETMX) , RHNX(NETMX)

C     Local variables.
      CHARACTER AIN*2
      CHARACTER*2 ATST(23)
      CHARACTER*6 HPOL(3) , PNET(6)
      CHARACTER*78 COM(5)
      LOGICAL ISOPEN
      INTEGER I , IB11 , IC11 , ID11 , IIEXK , IFLOW , IFRQ , IFTIMP , 
     &        IFTIMW , IGO , INC , IPED , IPTAG , IPTAGF , IPTAGT , 
     &        IPTAQ , IPTAQF , IPTAQT , IPTFLG , IPTFLQ , IRNGF , 
     &        ISAVE , ISEGNO , ITMP1 , ITMP2 , ITMP3 , ITMP4 , ITMP5 , 
     &        IX11 , J , JUMP , LAIN , MHZ , MPCNT , NFRQ , NPHI , 
     &        NPHIC , NTHI , NTHIC , NLOAD , NLODF , KCOM , IPSYM , 
     &        N1 , N2 , N , NP , M1 , M2 , M , MP , NVQD , NSANT , 
     &        NQDS , NCOUP , ICOUP , NEQ , NPEQ , NEQ2 , NONET , NTSOL , 
     &        NPRINT , MASYM , IPLP1 , IPLP2 , IPLP3 , IPLP4 , NTH , 
     &        NPH , IPD ,IAVP , INOR , IAX , IXTYP , NEAR , NFEH , NRX , 
     &        NRY , NRZ , JSNO , NSCON , NPCON , LOADSM
      REAL*8 CMAG , CVEL , DELFRQ , EPHA , EPHM , ETHA , ETHM , EXTIM , 
     &       FMHZ1 , FR , FR2 , PH , PHISS , RRKH , TA , THETIS , 
     &       TIM , TIM1 , TIM2 , TMP1 , TMP2 , TMP3 , TMP4 , TMP5 , 
     &       TMP6 , XPR1 , XPR2 , XPR3 , XPR4 , XPR5 , ZPNORM , EPSR , 
     &       SIG , SCRWLT , SCRWRT , FMHZ , WLAM , PIN , PNLS  , THETS , 
     &       PHIS ,DTH , DPH , RFLD , GNOR , CLT , CHT , EPSR2 , SIG2 , 
     &       XPR6 , PINR , PNLR , PLOSS , XNR , YNR , ZNR , DXNR , 
     &       DYNR , DZNR , OSIG , OEPSR , OFMHZ
      COMPLEX*16 FJ , ETH , EPH , CURI , EX , EY , EZ , EPSC , ZPED 

C     Common storage.
      INCLUDE 'matpar.inc'
      INCLUDE 'gnd.inc'
      INCLUDE 'ggrid.inc'
      INCLUDE 'gwav.inc'
 
C     Data initialisation.
      DATA ATST /'CE' , 'FR' , 'LD' , 'GN' , 'EX' , 'NT' , 'XQ' , 'NE' , 
     &           'GD' , 'RP' , 'CM' , 'NX' , 'EN' , 'TL' , 'PT' , 'KH' , 
     &           'NH' , 'PQ' , 'EK' , 'WG' , 'CP' , 'PL' , 'IG'/
      DATA HPOL /'LINEAR' , 'RIGHT' , 'LEFT'/
      DATA PNET /'      ' , '  ' , 'STRAIG' , 'HT' , 'CROSSE' , 'D'/
      DATA TA /0.0174532925199432957692369D0/
      DATA CVEL /299.792458D0/
      DATA OFMHZ /0.0D0/
      DATA OEPSR /0.0D0/
      DATA OSIG /0.0D0/
      DATA LOADSM /0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING NEC2D'
      ENDIF

      NRECNT = 1
      NRHCNT = 1 
      RDPCNT = 1

      CALL FILSUF(JOBNAM,FILBUF,0,SUFINP,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      OPEN (UNIT=CHINPT,FILE=FILBUF,STATUS='OLD',ERR=702)
      CALL FILSUF(JOBNAM,FILBUF,0,SUFRES,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      OPEN (UNIT=CHRSLT,FILE=FILBUF,STATUS='NEW',ERR=703)
      GOTO 705

  702 IFAIL=34
      RETURN

  703 IFAIL=49
      RETURN

  705 CONTINUE

      IF ( MULTIF.GT.0 ) THEN 
         CALL FILSUF(JOBNAM,FILBUF,0,SUFIMP,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         OPEN(UNIT=CHIMPD,STATUS='NEW',FILE=FILBUF)
      ENDIF

      CALL TIMER(EXTIM)
      FJ = (0.0D0,1.0D0)
C     IDF. This causes built in somnec to fail.
C      NXA(1) = 0
C     IDF
    1 KCOM = 0

      IFTIMW = 0
      IFTIMP = 0

    2 KCOM = KCOM + 1
      IF ( KCOM.GT.5 ) KCOM = 5
      READ (CHINPT,125) AIN , COM(KCOM)
      CALL UPCASE(AIN,AIN,LAIN,DEBUG)
      IF ( KCOM.GT.1 ) GOTO 3
      WRITE (CHRSLT,126)
      WRITE (CHRSLT,127)
      WRITE (CHRSLT,128)
    3 WRITE (CHRSLT,129) COM(KCOM)
      IF ( AIN.EQ.ATST(11) ) GOTO 2
      IF ( AIN.EQ.ATST(1) ) GOTO 4
      WRITE (CHRSLT,130)
      IFAIL=35
      RETURN
    4 CONTINUE
      DO 5 I = 1 , LD
         ZARRAY(I) = (0.0D0,0.0D0)
    5 CONTINUE
      MPCNT = 0
      IMAT = 0

C     Set up geometry data in subroutine DATAGN.

      CALL DATAGN(LD,SALP,IRESRV,CM,NLODF,ZARRAY,KCOM,COM,IP,EPSR,
     &            SIG,SCRWLT,SCRWRT,FMHZ,IPSYM,N1,N2,N,NP,M1,M2,M,MP,
     &            ICON1,ICON2,ITAG,ICONX,WLAM,X,Y,Z,SI,BI,ALP,BET,X2,Y2,
     &            Z2,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IPLP1,IPLP2,JMAX,
     &            NSMAXX,NPMAX,NSCON,NPCON,JCO,ISCON,IPCON,JOBNAM,
     &            FILBUF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IFLOW = 1
      IF ( IMAT.EQ.0 ) GOTO 326

C     Core allocation for arrays B, C, and D for N.G.F. solution.

      NEQ = N1 + 2*M1
      NEQ2 = N - N1 + 2*(M-M1) + NSCON + 2*NPCON
      CALL FBNGF(NEQ,NEQ2,IB11,IC11,ID11,IX11,IRESRV,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 6
  326 NEQ = N + 2*M
      NEQ2 = 0
      IB11 = 1
      IC11 = 1
      ID11 = 1
      IX11 = 1
      ICASX = 0
    6 NPEQ = NP + 2*MP
      WRITE (CHRSLT,135)

C     Default values for input parameters and flags.

      IPLP1 = 0
      IPLP2 = 0
      IPLP3 = 0
      IPLP4 = 0

      IGO = 1
      NFRQ = 1
      RRKH = 1.0D0
      IIEXK = 0
      IXTYP = 0
      NLOAD = 0
      NONET = 0
      NEAR = -1
      IPTFLG = -2
      IPTFLQ = -1
      IFAR = -1
      ZRATI = (1.0D0,0.0D0)
      IPED = 0
      IRNGF = 0
      NCOUP = 0
      ICOUP = 0
      IF ( ICASX.GT.0 ) GOTO 14
      FMHZ = CVEL
      NLODF = 0
      KSYMP = 1
      NRADL = 0
      IPERF = 0

C     Main input section - standard read statement - jumps to appro-
C     priate section for specific parameter set up.

   14 CALL READMN(2,AIN,ITMP1,ITMP2,ITMP3,ITMP4,TMP1,TMP2,TMP3,TMP4,
     &            TMP5,TMP6,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
C IDF ignore comments.
      IF ( AIN.EQ.ATST(23) ) GOTO 14
C End IDF.
      MPCNT = MPCNT + 1
      WRITE (CHRSLT,137) MPCNT , AIN , ITMP1 , ITMP2 , ITMP3 , ITMP4 , 
     &                   TMP1 , TMP2 , TMP3 , TMP4 , TMP5 , TMP6
      IF ( AIN.EQ.ATST(2) ) GOTO 16
      IF ( AIN.EQ.ATST(3) ) GOTO 17
      IF ( AIN.EQ.ATST(4) ) GOTO 21
      IF ( AIN.EQ.ATST(5) ) GOTO 24
      IF ( AIN.EQ.ATST(6) ) GOTO 28
      IF ( AIN.EQ.ATST(14) ) GOTO 28
      IF ( AIN.EQ.ATST(15) ) GOTO 31
      IF ( AIN.EQ.ATST(18) ) GOTO 319
      IF ( AIN.EQ.ATST(7) ) GOTO 37
      IF ( AIN.EQ.ATST(8) ) GOTO 32
      IF ( AIN.EQ.ATST(17) ) GOTO 208
      IF ( AIN.EQ.ATST(9) ) GOTO 34
      IF ( AIN.EQ.ATST(10) ) GOTO 36
      IF ( AIN.EQ.ATST(16) ) GOTO 305
      IF ( AIN.EQ.ATST(19) ) GOTO 320
      IF ( AIN.EQ.ATST(12) ) GOTO 1
      IF ( AIN.EQ.ATST(20) ) GOTO 322
      IF ( AIN.EQ.ATST(21) ) GOTO 304
      IF ( AIN.EQ.ATST(22) ) GOTO 330
      IF ( AIN.NE.ATST(13) ) GOTO 15
      CALL TIMER(TMP1)
      TMP1 = TMP1 - EXTIM
      WRITE (CHRSLT,201) TMP1

C     Normal end. Close open files.
      INQUIRE(UNIT=CHNREF,OPENED=ISOPEN)
      IF ( ISOPEN.EQV..TRUE. ) THEN
         CLOSE(CHNREF)
      ENDIF
      INQUIRE(UNIT=CHNRHF,OPENED=ISOPEN)
      IF ( ISOPEN.EQV..TRUE. ) THEN
         CLOSE(CHNRHF)
      ENDIF
      INQUIRE(UNIT=CHRDPT,OPENED=ISOPEN)
      IF ( ISOPEN.EQV..TRUE. ) THEN
         CLOSE(CHRDPT)
      ENDIF
      INQUIRE(UNIT=CHIMPD,OPENED=ISOPEN)
      IF ( ISOPEN.EQV..TRUE. ) THEN
         CLOSE(CHIMPD)
      ENDIF
      RETURN

   15 WRITE (CHRSLT,138)
      IFAIL=36
      RETURN

C     Frequency parameters.

   16 IFRQ = ITMP1
      IF ( ICASX.EQ.0 ) GOTO 8
      WRITE (CHRSLT,303) AIN
      IFAIL=37
      RETURN
    8 NFRQ = ITMP2
      IF ( NFRQ.EQ.0 ) NFRQ = 1
      FMHZ = TMP1
      DELFRQ = TMP2
      IF ( IPED.EQ.1 ) ZPNORM = 0.0D0
      IGO = 1
      IFLOW = 1
      GOTO 14

C     Matrix integration limit.

  305 RRKH = TMP1
      IF ( IGO.GT.2 ) IGO = 2
      IFLOW = 1
      GOTO 14

C     Extended thin wire kernel option.

  320 IIEXK = 1
      IF ( ITMP1.EQ.-1 ) IIEXK = 0
      IF ( IGO.GT.2 ) IGO = 2
      IFLOW = 1
      GOTO 14

C     Maximum coupling between antennas.

  304 IF ( IFLOW.NE.2 ) NCOUP = 0
      ICOUP = 0
      IFLOW = 2
      IF ( ITMP2.EQ.0 ) GOTO 14
      NCOUP = NCOUP + 1
      IF ( NCOUP.GT.MXCOUP ) GOTO 312
      NCTAG(NCOUP) = ITMP1
      NCSEG(NCOUP) = ITMP2
      IF ( ITMP4.EQ.0 ) GOTO 14
      NCOUP = NCOUP + 1
      IF ( NCOUP.GT.MXCOUP ) GOTO 312
      NCTAG(NCOUP) = ITMP3
      NCSEG(NCOUP) = ITMP4
      GOTO 14
  312 WRITE (CHRSLT,313)
      IFAIL=38
      RETURN

C     Loading parameters.

   17 IF ( IFLOW.EQ.3 ) GOTO 18
      NLOAD = 0
      IFLOW = 3
      IF ( IGO.GT.2 ) IGO = 2
      IF ( ITMP1.EQ.(-1) ) GOTO 14
   18 NLOAD = NLOAD + 1
      IF ( NLOAD.LE.LOADMX ) GOTO 19
      WRITE (CHRSLT,139)
      IFAIL=39
      RETURN
   19 LDTYP(NLOAD) = ITMP1
      LDTAG(NLOAD) = ITMP2
      IF ( ITMP4.EQ.0 ) ITMP4 = ITMP3
      LDTAGF(NLOAD) = ITMP3
      LDTAGT(NLOAD) = ITMP4
      IF ( ITMP4.GE.ITMP3 ) GOTO 20
      WRITE (CHRSLT,140) NLOAD , ITMP3 , ITMP4
      IFAIL=40
      RETURN
   20 ZLR(NLOAD) = TMP1
      ZLI(NLOAD) = TMP2
      ZLC(NLOAD) = TMP3
      GOTO 14

C     Ground parameters under the antenna.

   21 IFLOW = 4
      IF ( ICASX.EQ.0 ) GOTO 10
      WRITE (CHRSLT,303) AIN
      IFAIL=37
      RETURN
   10 IF ( IGO.GT.2 ) IGO = 2
      IF ( ITMP1.NE.(-1) ) GOTO 22
      KSYMP = 1
      NRADL = 0
      IPERF = 0
      GOTO 14
   22 IPERF = ITMP1
      NRADL = ITMP2
      KSYMP = 2
      EPSR = TMP1
      SIG = TMP2
      IF ( NRADL.EQ.0 ) GOTO 23
      IF ( IPERF.NE.2 ) GOTO 314
      WRITE (CHRSLT,390)
      IFAIL=41
      RETURN
  314 SCRWLT = TMP3
      SCRWRT = TMP4
      GOTO 14
   23 EPSR2 = TMP3
      SIG2 = TMP4
      CLT = TMP5
      CHT = TMP6
      GOTO 14

C     Excitation parameters.

   24 IF ( IFLOW.EQ.5 ) GOTO 25
      NSANT = 0
      NVQD = 0
      IPED = 0
      IFLOW = 5
      IF ( IGO.GT.3 ) IGO = 3
   25 MASYM = ITMP4/10
      IF ( ITMP1.GT.0 .AND. ITMP1.NE.5 ) GOTO 27
      IXTYP = ITMP1
      NTSOL = 0
      IF ( IXTYP.EQ.0 ) GOTO 205
      NVQD = NVQD + 1
      IF ( NVQD.GT.NSMAX ) GOTO 206
      IVQD(NVQD) = ISEGNO(ITMP2,ITMP3,LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      VQD(NVQD) = DCMPLX(TMP1,TMP2)
      IF ( CDABS(VQD(NVQD)).LT.1.D-20 ) VQD(NVQD) = (1.0D0,0.0D0)
      GOTO 207
  205 NSANT = NSANT + 1
      IF ( NSANT.LE.NSMAX ) GOTO 26
  206 WRITE (CHRSLT,141)
      IFAIL=42
      RETURN
   26 ISANT(NSANT) = ISEGNO(ITMP2,ITMP3,LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      VSANT(NSANT) = DCMPLX(TMP1,TMP2)
      IF ( CDABS(VSANT(NSANT)).LT.1.D-20 ) VSANT(NSANT) = (1.0D0,0.0D0)
  207 IPED = ITMP4 - MASYM*10
      ZPNORM = TMP3
      IF ( IPED.EQ.1 .AND. ZPNORM.GT.0 ) IPED = 2
      GOTO 14
   27 IF ( IXTYP.EQ.0 .OR. IXTYP.EQ.5 ) NTSOL = 0
      IXTYP = ITMP1
      NTHI = ITMP2
      NPHI = ITMP3
      XPR1 = TMP1
      XPR2 = TMP2
      XPR3 = TMP3
      XPR4 = TMP4
      XPR5 = TMP5
      XPR6 = TMP6
      NSANT = 0
      NVQD = 0
      THETIS = XPR1
      PHISS = XPR2
      GOTO 14

C     Network parameters.

   28 IF ( IFLOW.EQ.6 ) GOTO 29
      NONET = 0
      NTSOL = 0
      IFLOW = 6
      IF ( IGO.GT.3 ) IGO = 3
      IF ( ITMP2.EQ.(-1) ) GOTO 14
   29 NONET = NONET + 1
      IF ( NONET.LE.NETMX ) GOTO 30
      WRITE (CHRSLT,142)
      IFAIL=43
      RETURN
   30 NTYP(NONET) = 2
      IF ( AIN.EQ.ATST(6) ) NTYP(NONET) = 1
      ISEG1(NONET) = ISEGNO(ITMP1,ITMP2,LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      ISEG2(NONET) = ISEGNO(ITMP3,ITMP4,LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      X11R(NONET) = TMP1
      X11I(NONET) = TMP2
      X12R(NONET) = TMP3
      X12I(NONET) = TMP4
      X22R(NONET) = TMP5
      X22I(NONET) = TMP6
      IF ( NTYP(NONET).EQ.1 .OR. TMP1.GT.0. ) GOTO 14
      NTYP(NONET) = 3
      X11R(NONET) = -TMP1
      GOTO 14

C     Plot flags.

  330 IPLP1 = ITMP1
      IPLP2 = ITMP2
      IPLP3 = ITMP3
      IPLP4 = ITMP4

      GOTO 14

C     Print control for current.

   31 IPTFLG = ITMP1
      IPTAG = ITMP2
      IPTAGF = ITMP3
      IPTAGT = ITMP4
      IF ( ITMP3.EQ.0 .AND. IPTFLG.NE.-1 ) IPTFLG = -2
      IF ( ITMP4.EQ.0 ) IPTAGT = IPTAGF
      GOTO 14

C     Write control for charge.

  319 IPTFLQ = ITMP1
      IPTAQ = ITMP2
      IPTAQF = ITMP3
      IPTAQT = ITMP4
      IF ( ITMP3.EQ.0 .AND. IPTFLQ.NE.-1 ) IPTFLQ = -2
      IF ( ITMP4.EQ.0 ) IPTAQT = IPTAQF
      GOTO 14

C     Near field calculation parameters.

  208 NFEH = 1
      GOTO 209
   32 NFEH = 0
  209 IF ( IFLOW.NE.8 .OR. NFRQ.EQ.1 ) GOTO 33
      WRITE (CHRSLT,143)
   33 NEAR = ITMP1
      NRX = ITMP2
      NRY = ITMP3
      NRZ = ITMP4
      XNR = TMP1
      YNR = TMP2
      ZNR = TMP3
      DXNR = TMP4
      DYNR = TMP5
      DZNR = TMP6
      IFLOW = 8
      IF ( NFRQ.NE.1 ) GOTO 14
      GOTO (41,46,53,71,72) , IGO

C     Ground representation.

   34 EPSR2 = TMP1
      SIG2 = TMP2
      CLT = TMP3
      CHT = TMP4
      IFLOW = 9
      GOTO 14

C     Standard observation angle parameters.

   36 IFAR = ITMP1
      NTH = ITMP2
      NPH = ITMP3
      IF ( NTH.EQ.0 ) NTH = 1
      IF ( NPH.EQ.0 ) NPH = 1
      IPD = ITMP4/10
      IAVP = ITMP4 - IPD*10
      INOR = IPD/10
      IPD = IPD - INOR*10
      IAX = INOR/10
      INOR = INOR - IAX*10
      IF ( IAX.NE.0 ) IAX = 1
      IF ( IPD.NE.0 ) IPD = 1
      IF ( NTH.LT.2 .OR. NPH.LT.2 ) IAVP = 0
      IF ( IFAR.EQ.1 ) IAVP = 0
      THETS = TMP1
      PHIS = TMP2
      DTH = TMP3
      DPH = TMP4
      RFLD = TMP5
      GNOR = TMP6
      IFLOW = 10
      GOTO (41,46,53,71,78) , IGO

C     Write numerical green's function tape.

  322 IFLOW = 12
      IF ( ICASX.EQ.0 ) GOTO 301
      WRITE (CHRSLT,302)
      IFAIL=44
      RETURN
  301 IRNGF = IRESRV/2
      GOTO (41,46,52,52,52) , IGO

C     Execute card  -  calc. including radiated fields.

   37 IF ( IFLOW.EQ.10 .AND. ITMP1.EQ.0 ) GOTO 14
      IF ( NFRQ.EQ.1 .AND. ITMP1.EQ.0 .AND. IFLOW.GT.7 ) GOTO 14
      IF ( ITMP1.NE.0 ) GOTO 39
      IF ( IFLOW.GT.7 ) GOTO 38
      IFLOW = 7
      GOTO 40
   38 IFLOW = 11
      GOTO 40
   39 IFAR = 0
      RFLD = 0.0D0
      IPD = 0
      IAVP = 0
      INOR = 0
      IAX = 0
      NTH = 91
      NPH = 1
      THETS = 0.0D0
      PHIS = 0.0D0
      DTH = 1.0D0
      DPH = 0.0D0
      IF ( ITMP1.EQ.2 ) PHIS = 90.0D0
      IF ( ITMP1.NE.3 ) GOTO 40
      NPH = 2
      DPH = 90.0D0
   40 GOTO (41,46,53,71,78) , IGO

C     End of the main input section.

C     Beginning of the frequency do loop.

   41 MHZ = 1

      IF ( N.EQ.0 .OR. IFTIMW.EQ.1 ) GOTO 406
      IFTIMW = 1
      DO 445 I = 1 , N
         XTEMP(I) = X(I)
         YTEMP(I) = Y(I)
         ZTEMP(I) = Z(I)
         SITEMP(I) = SI(I)
         BITEMP(I) = BI(I)
  445 CONTINUE
  406 IF ( M.EQ.0 .OR. IFTIMP.EQ.1 ) GOTO 407
      IFTIMP = 1
      J = LD + 1
      DO 545 I = 1 , M
         J = J - 1
         XTEMP(J) = X(J)
         YTEMP(J) = Y(J)
         ZTEMP(J) = Z(J)
         BITEMP(J) = BI(J)
  545 CONTINUE
  407 CONTINUE
      FMHZ1 = FMHZ

C     Core allocation for primary interacton matrix.

      IF ( IMAT.EQ.0 ) THEN
         CALL FBLOCK(NPEQ,NEQ,IRESRV,IRNGF,IPSYM,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
      ENDIF
   42 IF ( MHZ.EQ.1 ) GOTO 44
      IF ( IFRQ.EQ.1 ) GOTO 43
C     FMHZ=FMHZ+DELFRQ

      FMHZ = FMHZ1 + (MHZ-1)*DELFRQ
      GOTO 44
   43 FMHZ = FMHZ*DELFRQ
   44 FR = FMHZ/CVEL

      WLAM = CVEL/FMHZ
      WRITE (CHRSLT,145) FMHZ , WLAM
      WRITE (CHRSLT,196) RRKH
      IF ( IIEXK.EQ.1 ) WRITE (CHRSLT,321)

C     Frequency scaling of geometric parameters.

      IF ( N.EQ.0 ) GOTO 306
      DO 45 I = 1 , N

         X(I) = XTEMP(I)*FR
         Y(I) = YTEMP(I)*FR
         Z(I) = ZTEMP(I)*FR
         SI(I) = SITEMP(I)*FR
         BI(I) = BITEMP(I)*FR
   45 CONTINUE

  306 IF ( M.EQ.0 ) GOTO 307
      FR2 = FR*FR
      J = LD + 1
      DO 245 I = 1 , M
         J = J - 1

         X(J) = XTEMP(J)*FR
         Y(J) = YTEMP(J)*FR
         Z(J) = ZTEMP(J)*FR
         BI(J) = BITEMP(J)*FR2
  245 CONTINUE

  307 IGO = 2

C     Structure segment loading.

   46 WRITE (CHRSLT,146)
      IF ( NLOAD.NE.0 ) THEN
         CALL LOAD(LDTYP,LDTAG,LDTAGF,LDTAGT,ZLR,ZLI,ZLC,LD,NLOAD,
     &        ZARRAY,N1,N2,N,NP,M1,ITAG,WLAM,SI,BI,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
      ENDIF
      IF ( NLOAD.EQ.0 .AND. NLODF.EQ.0 ) WRITE (CHRSLT,147)
      IF ( NLOAD.EQ.0 .AND. NLODF.NE.0 ) WRITE (CHRSLT,327)

C     Ground parameter.

      WRITE (CHRSLT,148)
      IF ( KSYMP.EQ.1 ) GOTO 49
      FRATI = (1.0D0,0.0D0)
      IF ( IPERF.EQ.1 ) GOTO 48
      IF ( SIG.LT.0. ) SIG = -SIG/(59.9584916D0*WLAM)
      EPSC = DCMPLX(EPSR,-SIG*WLAM*59.9584916D0)
      ZRATI = 1.0D0/CDSQRT(EPSC)
      U = ZRATI
      U2 = U*U
      IF ( NRADL.EQ.0 ) GOTO 47
      SCRWL = SCRWLT/WLAM
      SCRWR = SCRWRT/WLAM
      T1 = FJ*2367.066370312157358387101D0/DBLE(NRADL)
      T2 = SCRWR*DBLE(NRADL)
      WRITE (CHRSLT,170) NRADL , SCRWLT , SCRWRT
      WRITE (CHRSLT,149)
   47 IF ( IPERF.EQ.2 ) GOTO 328
      WRITE (CHRSLT,391)
      GOTO 329
#if defined (BUILTIN_SOMNEC)
  328 IF ( ( DABS( FMHZ - OFMHZ ).GT.1D-30 )
     & .OR.( DABS( EPSR - OEPSR ).GT.1D-30 )
     & .OR.( DABS( SIG - OSIG ).GT.1D-30 ) ) THEN
         IF(DEBUG.GT.0) THEN
            WRITE(*,9000) FMHZ,EPSR,SIG
         ENDIF
         CALL SOMNEC(IFAIL,0,DEBUG,0,EPSR,SIG,FMHZ,'NOOUT')
         OFMHZ = FMHZ
         OEPSR = EPSR
         OSIG = SIG
         IF(IFAIL.NE.0) RETURN
      ENDIF
#else
 328  IF (LOADSM.EQ.0) THEN 
         CALL FILSUF(JOBNAM,FILBUF,0,SUFSOM,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         OPEN (UNIT=CHSOMN,FILE=FILBUF,STATUS='OLD',
     &         FORM='UNFORMATTED',ERR=800)
         GOTO 801
  800    WRITE (CHRSLT,900)
         IFAIL=45
         RETURN
  801    READ (CHSOMN) AR1 , AR2 , AR3 , EPSCF , DXA , DYA , XSA , YSA , 
     &             NXA , NYA
         LOADSM=1
      ENDIF
#endif
      FRATI = (EPSC-1.0D0)/(EPSC+1.0D0)
      IF ( CDABS((EPSCF-EPSC)/EPSC).LT.1.D-3 ) GOTO 400
      WRITE (CHRSLT,393) EPSCF , EPSC
      IFAIL=46
      RETURN
  400 WRITE (CHRSLT,392)
  329 WRITE (CHRSLT,150) EPSR , SIG , EPSC
      GOTO 50
   48 WRITE (CHRSLT,151)
      GOTO 50
   49 WRITE (CHRSLT,152)
   50 CONTINUE

C     Fill and factor primary interaction matrix.

      CALL TIMER(TIM1)
      IF ( ICASX.NE.0 ) GOTO 324
      CALL CMSET(NEQ,CM,RRKH,IIEXK,LD,SALP,NLOAD,ZARRAY,D,N1,N,NP,M1,M,
     &           MP,ICON1,ICON2,ICONX,X,Y,Z,SI,BI,CAB,SAB,T1X,T1Y,T1Z,
     &           T2X,T2Y,T2Z,JMAX,JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL TIMER(TIM2)
      TIM = TIM2 - TIM1
      CALL FACTRS(NPEQ,NEQ,CM,IP,IX,11,12,13,14,LD,D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      GOTO 323

C     N.G.F. - Fill B, C, and D and factor D-C(INV(A)B).

  324 IF ( NEQ2.EQ.0 ) GOTO 333

      CALL CMNGF(CM(IB11),CM(IC11),CM(ID11),NPBX,NEQ,NEQ2,RRKH,IIEXK,
     &           LD,SALP,NLOAD,NLODF,ZARRAY,N1,N2,N,NP,M1,M2,M,MP,ICON1,
     &           ICON2,ICONX,X,Y,Z,SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,
     &           T2Z,JMAX,NSMAXX,NPMAX,JSNO,NSCON,NPCON,JCO,ISCON,IPCON,
     &           AX,BX,CX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL TIMER(TIM2)
      TIM = TIM2 - TIM1
      CALL FACGF(CM,CM(IB11),CM(IC11),CM(ID11),CM(IX11),IP,IX,NP,N1,MP,
     &           M1,NEQ,NEQ2,LD,D,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
  323 CALL TIMER(TIM1)
      TIM2 = TIM1 - TIM2
      WRITE (CHRSLT,153) TIM , TIM2
  333 IGO = 3
      NTSOL = 0
      IF ( IFLOW.NE.12 ) GOTO 53

C     Write N.G.F. file.

   52 CALL GFOUT(LD,SALP,IRESRV,CM,NLOAD,NLODF,ZARRAY,KCOM,COM,IP,
     &           EPSR,SIG,SCRWLT,SCRWRT,FMHZ,IPSYM,N1,N2,N,NP,M1,M2,M,
     &           MP,ICON1,ICON2,ITAG,ICONX,WLAM,X,Y,Z,SI,BI,ALP,BET,T2X,
     &           T2Y,T2Z,JOBNAM,FILBUF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN

      GOTO 14

C     Excitation set up (right hand side, -Einc.).

   53 NTHIC = 1
      NPHIC = 1
      INC = 1
      NPRINT = 0
   54 IF ( IXTYP.EQ.0 .OR. IXTYP.EQ.5 ) GOTO 56
      IF ( IPTFLG.LE.0 .OR. IXTYP.EQ.4 ) WRITE (CHRSLT,154)
      TMP5 = TA*XPR5
      TMP4 = TA*XPR4
      IF ( IXTYP.NE.4 ) GOTO 55
      TMP1 = XPR1/WLAM
      TMP2 = XPR2/WLAM
      TMP3 = XPR3/WLAM
      TMP6 = XPR6/(WLAM*WLAM)
      WRITE (CHRSLT,156) XPR1 , XPR2 , XPR3 , XPR4 , XPR5 , XPR6
      GOTO 56
   55 TMP1 = TA*XPR1
      TMP2 = TA*XPR2
      TMP3 = TA*XPR3
      TMP6 = XPR6
      IF ( IPTFLG.LE.0 ) WRITE (CHRSLT,155) XPR1 , XPR2 , XPR3 , 
     &                                      HPOL(IXTYP) , XPR6
   56 CALL ETMNS(TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,IXTYP,CUR,LD,SALP,
     &           NLOAD,NLODF,ZARRAY,N,M,ICON1,ICON2,WLAM,X,Y,Z,SI,BI,
     &           CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NVQD,NSANT,NQDS,
     &           IVQD,ISANT,IQDS,VQD,VSANT,VQDS,JMAX,JSNO,JCO,AX,BX,CX,
     &           IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
C     Matrix solving (NETWK calls solves).

      IF ( NONET.EQ.0 .OR. INC.GT.1 ) GOTO 60
      WRITE (CHRSLT,158)
      ITMP3 = 0
      ITMP1 = NTYP(1)
      DO 59 I = 1 , 2
         IF ( ITMP1.EQ.3 ) ITMP1 = 2
         IF ( ITMP1.EQ.2 ) WRITE (CHRSLT,159)
         IF ( ITMP1.EQ.1 ) WRITE (CHRSLT,160)
         DO 58 J = 1 , NONET
            ITMP2 = NTYP(J)
            IF ( (ITMP2/ITMP1).EQ.1 ) GOTO 57
            ITMP3 = ITMP2
            GOTO 58
   57       ITMP4 = ISEG1(J)
            ITMP5 = ISEG2(J)
            IF ( ITMP2.GE.2 .AND. X11I(J).LE.0. ) X11I(J)
     &           = WLAM*DSQRT((X(ITMP5)-X(ITMP4))**2+(Y(ITMP5)-Y(ITMP4))
     &           **2+(Z(ITMP5)-Z(ITMP4))**2)
            WRITE (CHRSLT,157) ITAG(ITMP4) , ITMP4 , ITAG(ITMP5) , 
     &                         ITMP5 , X11R(J) , X11I(J) , X12R(J) , 
     &                         X12I(J) , X22R(J) , X22I(J) , 
     &                         PNET(2*ITMP2-1) , PNET(2*ITMP2)
   58    CONTINUE
         IF ( ITMP3.EQ.0 ) GOTO 60
         ITMP1 = ITMP3
   59 CONTINUE
   60 CONTINUE
      IF ( INC.GT.1 .AND. IPTFLG.GT.0 ) NPRINT = 1
      CALL NETWK(CM,CM(IB11),CM(IC11),CM(ID11),IP,CUR,FMHZ,LD,AIR,
     &           AII,BIR,BII,CIR,CII,D,N1,N,NP,M1,M,MP,ICON1,ICON2,
     &           ITAG,WLAM,SI,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NVQD,
     &           NSANT,NQDS,IVQD,ISANT,IQDS,VQD,VSANT,VQDS,NETMX,NEQ,
     &           NEQ2,NONET,NTSOL,NPRINT,MASYM,ISEG1,ISEG2,NTYP,PIN,
     &           PNLS,X11R,X11I,X12R,X12I,X22R,X22I,ZPED,JMAX,NSMAXX,
     &           JSNO,NSCON,NPCON,JCO,ISCON,AX,BX,CX,MULTIF,IPNT,NTEQA, 
     &           NTSCA,CMN,RHNT,RHS,VSRC,RHNX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      NTSOL = 1
      IF ( IPED.EQ.0 ) GOTO 61
      ITMP1 = MHZ + 4*(MHZ-1)
      IF ( ITMP1.GT.(NORMF-3) ) GOTO 61
      FNORM(ITMP1) = DBLE(ZPED)
      FNORM(ITMP1+1) = DIMAG(ZPED)
      FNORM(ITMP1+2) = CDABS(ZPED)
      FNORM(ITMP1+3) = CANG(ZPED)
      IF ( IPED.EQ.2 ) GOTO 61
      IF ( FNORM(ITMP1+2).GT.ZPNORM ) ZPNORM = FNORM(ITMP1+2)
   61 CONTINUE

C     Printing structure currents.

      IF ( N.EQ.0 ) GOTO 308
      IF ( IPTFLG.EQ.(-1) ) GOTO 63
      IF ( IPTFLG.GT.0 ) GOTO 62
      WRITE (CHRSLT,161)
      WRITE (CHRSLT,162)
      GOTO 63
   62 IF ( IPTFLG.EQ.3 .OR. INC.GT.1 ) GOTO 63
      WRITE (CHRSLT,163) XPR3 , HPOL(IXTYP) , XPR6
   63 PLOSS = 0.0D0
      ITMP1 = 0
      JUMP = IPTFLG + 1
      DO 69 I = 1 , N
         CURI = CUR(I)*WLAM
         CMAG = CDABS(CURI)
         PH = CANG(CURI)
         IF ( NLOAD.EQ.0 .AND. NLODF.EQ.0 ) GOTO 64
         IF ( DABS(DBLE(ZARRAY(I))).LT.1.D-20 ) GOTO 64
         PLOSS = PLOSS + 0.5D0*CMAG*CMAG*DBLE(ZARRAY(I))*SI(I)
   64    IF ( JUMP ) 68 , 69 , 65
   65    IF ( IPTAG.EQ.0 ) GOTO 66
         IF ( ITAG(I).NE.IPTAG ) GOTO 69
   66    ITMP1 = ITMP1 + 1
         IF ( ITMP1.LT.IPTAGF .OR. ITMP1.GT.IPTAGT ) GOTO 69
         IF ( IPTFLG.EQ.0 ) GOTO 68
         IF ( IPTFLG.LT.2 .OR. INC.GT.NORMF ) GOTO 67
         FNORM(INC) = CMAG
         ISAVE = I
   67    IF ( IPTFLG.NE.3 ) WRITE (CHRSLT,164) XPR1 , XPR2 , CMAG , 
     &                                         PH , I
         GOTO 69
   68    WRITE (CHRSLT,165) I , ITAG(I) , X(I) , Y(I) , Z(I) , SI(I) , 
     &                      CURI , CMAG , PH

         IF ( IPLP1.NE.1 ) GOTO 69
         IF ( IPLP2.EQ.1 ) WRITE (CHRPAT,*) CURI
         IF ( IPLP2.EQ.2 ) WRITE (CHRPAT,*) CMAG , PH

   69 CONTINUE
      IF ( IPTFLQ.EQ.(-1) ) GOTO 308
      WRITE (CHRSLT,315)
      ITMP1 = 0
      FR = 1.D-6/FMHZ
      DO 316 I = 1 , N
         IF ( IPTFLQ.EQ.(-2) ) GOTO 318
         IF ( IPTAQ.EQ.0 ) GOTO 317
         IF ( ITAG(I).NE.IPTAQ ) GOTO 316
  317    ITMP1 = ITMP1 + 1
         IF ( ITMP1.LT.IPTAQF .OR. ITMP1.GT.IPTAQT ) GOTO 316
  318    CURI = FR*DCMPLX(-BII(I),BIR(I))
         CMAG = CDABS(CURI)
         PH = CANG(CURI)
         WRITE (CHRSLT,165) I , ITAG(I) , X(I) , Y(I) , Z(I) , SI(I) , 
     &                      CURI , CMAG , PH
  316 CONTINUE
  308 IF ( M.EQ.0 ) GOTO 310
      WRITE (CHRSLT,197)
      J = N - 2
      ITMP1 = LD + 1
      DO 309 I = 1 , M
         J = J + 3
         ITMP1 = ITMP1 - 1
         EX = CUR(J)
         EY = CUR(J+1)
         EZ = CUR(J+2)
         ETH = EX*T1X(ITMP1) + EY*T1Y(ITMP1) + EZ*T1Z(ITMP1)
         EPH = EX*T2X(ITMP1) + EY*T2Y(ITMP1) + EZ*T2Z(ITMP1)
         ETHM = CDABS(ETH)
         ETHA = CANG(ETH)
         EPHM = CDABS(EPH)
         EPHA = CANG(EPH)

         WRITE (CHRSLT,198) I , X(ITMP1) , Y(ITMP1) , Z(ITMP1) , ETHM , 
     &                 ETHA , EPHM , EPHA , EX , EY , EZ
         IF ( IPLP1.NE.1 ) GOTO 309
         IF ( IPLP3.EQ.1 ) WRITE (CHRPAT,*) EX
         IF ( IPLP3.EQ.2 ) WRITE (CHRPAT,*) EY
         IF ( IPLP3.EQ.3 ) WRITE (CHRPAT,*) EZ
         IF ( IPLP3.EQ.4 ) WRITE (CHRPAT,*) EX , EY , EZ
  309 CONTINUE

  310 IF ( IXTYP.NE.0 .AND. IXTYP.NE.5 ) GOTO 70
      TMP1 = PIN - PNLS - PLOSS
      TMP2 = 100.0D0*TMP1/PIN
      WRITE (CHRSLT,166) PIN , TMP1 , PLOSS , PNLS , TMP2
   70 CONTINUE
      IGO = 4
      IF ( NCOUP.GT.0 ) CALL COUPLE(CUR,LD,N,ITAG,WLAM,NSMAX,NVQD,
     &                              NSANT,ISANT,VSANT,MXCOUP,MXCOU2,
     &                              NCOUP,ICOUP,NCTAG,NCSEG,Y11A,Y12A,
     &                              IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IF ( IFLOW.NE.7 ) GOTO 71
      IF ( IXTYP.GT.0 .AND. IXTYP.LT.4 ) GOTO 113
      IF ( NFRQ.NE.1 ) GOTO 120
      WRITE (CHRSLT,135)
      GOTO 14
   71 IGO = 5

C     Near field calculation.

   72 IF ( NEAR.EQ.(-1) ) GOTO 78
      CALL NFPAT(LD,SALP,AIR,AII,BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,
     &           WLAM,X,Y,Z,SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IPLP1,
     &           IPLP2,IPLP3,IPLP4,NEAR,NFEH,NRX,NRY,NRZ,XNR,YNR,ZNR,
     &           DXNR,DYNR,DZNR,JOBNAM,FILBUF,MULTIF,NRECNT,NRHCNT,
     &           IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IF ( MHZ.EQ.NFRQ ) NEAR = -1
      IF ( NFRQ.NE.1 ) GOTO 78
      WRITE (CHRSLT,135)
      GOTO 14

C     Standard far field calculation.

   78 IF ( IFAR.EQ.-1 ) GOTO 113
      PINR = PIN
      PNLR = PNLS
      CALL RDPAT(LD,SALP,AIR,AII,BIR,BII,CIR,CII,CUR,SCRWLT,
     &           SCRWRT,N,M,WLAM,X,Y,Z,SI,GAIN,CAB,SAB,XS,YS,ZS,Z2S,
     &           IPLP1,IPLP2,IPLP3,IPLP4,NTH,NPH,IPD,IAVP,INOR,IAX,
     &           IXTYP,THETS,PHIS,DTH,DPH,RFLD,GNOR,CLT,CHT,EPSR2,SIG2,
     &           XPR6,PINR,PNLR,PLOSS,JOBNAM,FILBUF,MULTIF,RDPCNT,
     &           IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
  113 IF ( IXTYP.EQ.0 .OR. IXTYP.GE.4 ) GOTO 119
      NTHIC = NTHIC + 1
      INC = INC + 1
      XPR1 = XPR1 + XPR4
      IF ( NTHIC.LE.NTHI ) GOTO 54
      NTHIC = 1
      XPR1 = THETIS
      XPR2 = XPR2 + XPR5
      NPHIC = NPHIC + 1
      IF ( NPHIC.LE.NPHI ) GOTO 54
      NPHIC = 1
      XPR2 = PHISS
      IF ( IPTFLG.LT.2 ) GOTO 119
C     Normalized receiving pattern printed.
      ITMP1 = NTHI*NPHI
      IF ( ITMP1.LE.NORMF ) GOTO 114
      ITMP1 = NORMF
      WRITE (CHRSLT,181)
  114 TMP1 = FNORM(1)
      DO 115 J = 2 , ITMP1
         IF ( FNORM(J).GT.TMP1 ) TMP1 = FNORM(J)
  115 CONTINUE
      WRITE (CHRSLT,182) TMP1 , XPR3 , HPOL(IXTYP) , XPR6 , ISAVE
      DO 118 J = 1 , NPHI
         ITMP2 = NTHI*(J-1)
         DO 116 I = 1 , NTHI
            ITMP3 = I + ITMP2
            IF ( ITMP3.GT.ITMP1 ) GOTO 117
            TMP2 = FNORM(ITMP3)/TMP1
            TMP3 = DB20(TMP2)
            WRITE (CHRSLT,183) XPR1 , XPR2 , TMP3 , TMP2
            XPR1 = XPR1 + XPR4
  116    CONTINUE
  117    XPR1 = THETIS
         XPR2 = XPR2 + XPR5
  118 CONTINUE
      XPR2 = PHISS
  119 IF ( MHZ.EQ.NFRQ ) IFAR = -1
      IF ( NFRQ.NE.1 ) GOTO 120
      WRITE (CHRSLT,135)
      GOTO 14
  120 MHZ = MHZ + 1
      IF ( MHZ.LE.NFRQ ) GOTO 42
      IF ( IPED.EQ.0 ) GOTO 123
      IF ( NVQD.LT.1 ) GOTO 199
      WRITE (CHRSLT,184) IVQD(NVQD) , ZPNORM
      GOTO 204
  199 WRITE (CHRSLT,184) ISANT(NSANT) , ZPNORM
  204 ITMP1 = NFRQ
      IF ( ITMP1.LE.(NORMF/4) ) GOTO 121
      ITMP1 = NORMF/4
      WRITE (CHRSLT,185)
  121 IF ( IFRQ.EQ.0 ) TMP1 = FMHZ - (NFRQ-1)*DELFRQ
      IF ( IFRQ.EQ.1 ) TMP1 = FMHZ/(DELFRQ**(NFRQ-1))
      DO 122 I = 1 , ITMP1
         ITMP2 = I + 4*(I-1)
         TMP2 = FNORM(ITMP2)/ZPNORM
         TMP3 = FNORM(ITMP2+1)/ZPNORM
         TMP4 = FNORM(ITMP2+2)/ZPNORM
         TMP5 = FNORM(ITMP2+3)
         WRITE (CHRSLT,186) TMP1 , FNORM(ITMP2) , FNORM(ITMP2+1) , 
     &                 FNORM(ITMP2+2) , FNORM(ITMP2+3) , TMP2 , TMP3 , 
     &                 TMP4 , TMP5
         IF ( IFRQ.EQ.0 ) TMP1 = TMP1 + DELFRQ
         IF ( IFRQ.EQ.1 ) TMP1 = TMP1*DELFRQ
  122 CONTINUE
      WRITE (CHRSLT,135)
  123 CONTINUE
      NFRQ = 1
      MHZ = 1
      GOTO 14
 
  125 FORMAT (A2,A)
  126 FORMAT ('1')
  127 FORMAT (///,33X,'*********************************************',
     &        //,36X,'NUMERICAL ELECTROMAGNETICS CODE (NEC-2D)',//,33X,
     &        '*********************************************')
  128 FORMAT (////,37X,'- - - - COMMENTS - - - -',//)
  129 FORMAT (25X,A)
  130 FORMAT (///,10X,'INCORRECT LABEL FOR A COMMENT CARD')
  135 FORMAT (/////)
C 136 FORMAT (A2,I3,3I5,6E10.3)
  137 FORMAT (1X,'***** DATA CARD NO.',I3,3X,A2,1X,I3,3(1X,I5),
     &        6(1X,1P,E12.5))
  138 FORMAT (///,10X,'FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION')
  139 FORMAT (///,10X,'NUMBER OF LOADING CARDS EXCEEDS STORAGE ALLOTTED'
     &        )
  140 FORMAT (///,10X,'DATA FAULT ON LOADING CARD NO.=',I5,5X,
     &        'ITAG STEP1=',I5,'  IS GREATER THAN ITAG STEP2=',I5)
  141 FORMAT (///,10X,
     &        'NUMBER OF EXCITATION CARDS EXCEEDS STORAGE ALLOTTED')
  142 FORMAT (///,10X,'NUMBER OF NETWORK CARDS EXCEEDS STORAGE ALLOTTED'
     &        )
  143 FORMAT (///,10X,
     &'WHEN MULTIPLE FREQUENCIES ARE REQUESTED, ONLY ONE NEAR FIELD CARD
     & CAN BE USED -',/,10X,'LAST CARD READ IS USED')
  145 FORMAT (////,33X,'- - - - - - FREQUENCY - - - - - -',//,36X,
     &        'FREQUENCY=',1P,E11.4,' MHZ',/,36X,'WAVELENGTH=',E11.4,
     &        ' METERS')
  146 FORMAT (///,30X,' - - - STRUCTURE IMPEDANCE LOADING - - -')
  147 FORMAT (/,35X,'THIS STRUCTURE IS NOT LOADED')
  148 FORMAT (///,34X,'- - - ANTENNA ENVIRONMENT - - -',/)
  149 FORMAT (40X,'MEDIUM UNDER SCREEN -')
  150 FORMAT (40X,'RELATIVE DIELECTRIC CONST.=',F7.3,/,40X,
     &        'CONDUCTIVITY=',1P,E10.3,' MHOS/METER',/,40X,
     &        'COMPLEX DIELECTRIC CONSTANT=',2E12.5)
  151 FORMAT (42X,'PERFECT GROUND')
  152 FORMAT (44X,'FREE SPACE')
  153 FORMAT (///,32X,'- - - MATRIX TIMING - - -',//,24X,'FILL=',F9.3,
     &        ' SEC.,  FACTOR=',F9.3,' SEC.')
  154 FORMAT (///,40X,'- - - EXCITATION - - -')
  155 FORMAT (/,4X,'PLANE WAVE',4X,'THETA=',F7.2,' DEG,  PHI=',F7.2,
     &        ' DEG,  ETA=',F7.2,' DEG,  TYPE -',A6,'=  AXIAL RATIO=',
     &        F6.3)
  156 FORMAT (/,31X,'POSITION (METERS)',14X,'ORIENTATION (DEG)='/,28X,
     &        'X',12X,'Y',12X,'Z',10X,'ALPHA',5X,'BETA',4X,
     &        'DIPOLE MOMENT',//,4X,'CURRENT SOURCE',1X,3(3X,F10.5),1X,
     &        2(3X,F7.2),4X,F8.3)
  157 FORMAT (4X,4(I5,1X),1P,6(3X,E11.4),3X,A6,A2)
  158 FORMAT (///,44X,'- - - NETWORK DATA - - -')
  159 FORMAT (/,6X,'- FROM -    - TO -',11X,'TRANSMISSION LINE',15X,
     &        '-  -  SHUNT ADMITTANCES (MHOS)  -  -',14X,'LINE',/,6X,
     &        'TAG  SEG.   TAG  SEG.',6X,'IMPEDANCE',6X,'LENGTH',12X,
     &        '- END ONE -',17X,'- END TWO -',12X,'TYPE',/,6X,
     &        'NO.   NO.   NO.   NO.',9X,'OHMS',8X,'METERS',9X,'REAL',
     &        10X,'IMAG.',9X,'REAL',10X,'IMAG.')
  160 FORMAT (/,6X,'- FROM -',4X,'- TO -',26X,
     &        '-  -  ADMITTANCE MATRIX ELEMENTS (MHOS)  -  -',/,6X,
     &        'TAG  SEG.   TAG  SEG.',13X,'(ONE,ONE)',19X,'(ONE,TWO)',
     &        19X,'(TWO,TWO)',/,6X,'NO.   NO.   NO.   NO.',8X,'REAL',
     &        10X,'IMAG.',9X,'REAL',10X,'IMAG.',9X,'REAL',10X,'IMAG.')
  161 FORMAT (///,29X,'- - - CURRENTS AND LOCATION - - -',//,33X,
     &        'DISTANCES IN WAVELENGTHS')
  162 FORMAT (//,2X,'SEG.',2X,'TAG',4X,'COORD. OF SEG. CENTER',5X,
     &        'SEG.',12X,'- - - CURRENT (AMPS) - - -',/,2X,'NO.',3X,
     &        'NO.',5X,'X',8X,'Y',8X,'Z',6X,'LENGTH',5X,'REAL',8X,
     &        'IMAG.',7X,'MAG.',8X,'PHASE')
  163 FORMAT (///,33X,'- - - RECEIVING PATTERN PARAMETERS - - -',/,43X,
     &        'ETA=',F7.2,' DEGREES',/,43X,'TYPE -',A6,/,43X,
     &        'AXIAL RATIO=',F6.3,//,11X,'THETA',6X,'PHI',10X,
     &        '-  CURRENT  -',9X,'SEG',/,11X,'(DEG)',5X,'(DEG)',7X,
     &        'MAGNITUDE',4X,'PHASE',6X,'NO.',/)
  164 FORMAT (10X,2(F7.2,3X),1X,1P,E11.4,3X,0P,F7.2,4X,I5)
  165 FORMAT (1X,2I5,3F9.4,F9.5,1X,1P,3E12.4,0P,F9.3)
  166 FORMAT (///,40X,'- - - POWER BUDGET - - -',//,43X,
     &        'INPUT POWER   =',1P,E11.4,' WATTS',/,43X,
     &        'RADIATED POWER=',E11.4,' WATTS',/,43X,'STRUCTURE LOSS=',
     &        E11.4,' WATTS',/,43X,'NETWORK LOSS  =',E11.4,' WATTS',/,
     &        43X,'EFFICIENCY    =',0P,F7.2,' PERCENT')
  170 FORMAT (40X,'RADIAL WIRE GROUND SCREEN',/,40X,I5,' WIRES',/,40X,
     &        'WIRE LENGTH=',F8.2,' METERS',/,40X,'WIRE RADIUS=',1P,
     &        E10.3,' METERS')
  181 FORMAT (///,4X,
     &        'RECEIVING PATTERN STORAGE TOO SMALL,ARRAY TRUNCATED')
  182 FORMAT (///,32X,'- - - NORMALIZED RECEIVING PATTERN - - -',/,41X,
     &        'NORMALIZATION FACTOR=',1P,E11.4,/,41X,'ETA=',0P,F7.2,
     &        ' DEGREES',/,41X,'TYPE -',A6,/,41X,'AXIAL RATIO=',F6.3,/,
     &        41X,'SEGMENT NO.=',I5,//,21X,'THETA',6X,'PHI',9X,
     &        '-  PATTERN  -',/,21X,'(DEG)',5X,'(DEG)',8X,'DB',8X,
     &        'MAGNITUDE',/)
  183 FORMAT (20X,2(F7.2,3X),1X,F7.2,4X,1P,E11.4)
  184 FORMAT (///,36X,'- - - INPUT IMPEDANCE DATA - - -',/,45X,
     &        'SOURCE SEGMENT NO.',I4,/,45X,'NORMALIZATION FACTOR=',1P,
     &        E12.5,//,7X,'FREQ.',13X,
     &        '-  -  UNNORMALIZED IMPEDANCE  -  -',21X,
     &        '-  -  NORMALIZED IMPEDANCE  -  -',/,19X,'RESISTANCE',4X,
     &        'REACTANCE',6X,'MAGNITUDE',4X,'PHASE',7X,'RESISTANCE',4X,
     &        'REACTANCE',6X,'MAGNITUDE',4X,'PHASE',/,8X,'MHZ',11X,
     &        'OHMS',10X,'OHMS',11X,'OHMS',5X,'DEGREES',47X,'DEGREES',/)
  185 FORMAT (///,4X,
     &  'STORAGE FOR IMPEDANCE NORMALIZATION TOO SMALL, ARRAY TRUNCATED'
     &  )
  186 FORMAT (3X,F9.3,2X,1P,2(2X,E12.5),3X,E12.5,2X,0P,F7.2,2X,1P,
     &        2(2X,E12.5),3X,E12.5,2X,0P,F7.2)
  196 FORMAT (////,20X,
     &        'APPROXIMATE INTEGRATION EMPLOYED FOR SEGMENTS MORE THAN',
     &        F8.3,' WAVELENGTHS APART')
  197 FORMAT (////,41X,'- - - - SURFACE PATCH CURRENTS - - - -',//,50X,
     &        'DISTANCE IN WAVELENGTHS',/,50X,'CURRENT IN AMPS/METER',
     &        //,28X,'- - SURFACE COMPONENTS - -',19X,
     &        '- - - RECTANGULAR COMPONENTS - - -',/,6X,'PATCH CENTER',
     &        6X,'TANGENT VECTOR 1',3X,'TANGENT VECTOR 2',11X,'X',19X,
     &        'Y',19X,'Z',/,5X,'X',6X,'Y',6X,'Z',5X,'MAG.',7X,'PHASE',
     &        3X,'MAG.',7X,'PHASE',3(4X,'REAL',6X,'IMAG. '))
  198 FORMAT (1X,I4,/,1X,3F7.3,2(1P,E11.4,0P,F8.2),1P,6E10.2)
  201 FORMAT (/,' RUN TIME =',F10.3)
  315 FORMAT (///,34X,'- - - CHARGE DENSITIES - - -',//,36X,
     &        'DISTANCES IN WAVELENGTHS',///,2X,'SEG.',2X,'TAG',4X,
     &        'COORD. OF SEG. CENTER',5X,'SEG.',10X,
     &        'CHARGE DENSITY (COULOMBS/METER)',/,2X,'NO.',3X,'NO.',5X,
     &        'X',8X,'Y',8X,'Z',6X,'LENGTH',5X,'REAL',8X,'IMAG.',7X,
     &        'MAG.',8X,'PHASE')
  321 FORMAT (/,20X,'THE EXTENDED THIN WIRE KERNEL WILL BE USED')
  303 FORMAT (/,' ERROR - ',A2,' CARD IS NOT ALLOWED WITH N.G.F.')
  327 FORMAT (/,35X,' LOADING ONLY IN N.G.F. SECTION')
  302 FORMAT (' ERROR - N.G.F. IN USE.  CANNOT WRITE NEW N.G.F.')
  313 FORMAT (/,
     &  ' NUMBER OF SEGMENTS IN COUPLING CALCULATION (CP) EXCEEDS LIMIT'
     &  )
  390 FORMAT (
     &' RADIAL WIRE G. S. APPROXIMATION MAY NOT BE USED WITH SOMMERFELD 
     &GROUND OPTION')
  391 FORMAT (40X,'FINITE GROUND.  REFLECTION COEFFICIENT APPROXIMATION'
     &        )
  392 FORMAT (40X,'FINITE GROUND.  SOMMERFELD SOLUTION')
  393 FORMAT (/,' ERROR IN GROUND PARAMETERS -',/,
     &        ' COMPLEX DIELECTRIC CONSTANT FROM FILE IS',1P,2E12.5,/,
     &        32X,'REQUESTED',2E12.5)
  900 FORMAT (' ERROR OPENING SOMMERFELD GROUND FILE - SOM2D.NEC')
 9000 FORMAT (' SOMNEC CALLED WITH F= ',E12.6,' ER= ',E12.6,' SIG= ',
     &        E12.6)

      END
