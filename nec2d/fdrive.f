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

      PROGRAM FDRIVE

C*--------------------------------------------------------------------**
C*                                                                    **
C* Fortran driver fron NEC2D using static allocation with limits in   **
C* limits.inc.                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ERROR   FILSUF  NEC2D                                  **
C* called by   ** NOTHING **                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL NEC2D , FILSUF

C     Parameter definitions.
      INCLUDE 'limits.inc'
      INCLUDE 'nec2d.inc'

C     Local variables.
      CHARACTER*1024 INFILE , OUTFIL , JOBNAM , REST
      CHARACTER*4    DOTSUF
      CHARACTER*1044 FILBUF
      INTEGER MULTIF , IFAIL , ISUF , DEBUG
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
     &           Y11A(MXCOUP) , Y12A(MXCOU2) , CMN(NETMX,NETMX) , 
     &           RHNT(NETMX) , RHS(3*LD) , VSRC(NETMX) , RHNX(NETMX)

C     Equivalent local storage.
      EQUIVALENCE (X2,SI) , (Y2,ALP) , (Z2,BET) , (CAB,ALP) , 
     &            (SAB,BET) , (T1X,SI) , (T1Y,ALP) , (T1Z,BET) ,
     &            (T2X,ICON1) , (T2Y,ICON2) , (T2Z,ITAG) ,
     &            (XS,X) , (YS,Y) , (ZS,Z) , (Z2S,BI)


      MULTIF = 2
      DEBUG  = 0

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FDRIVE'
      ENDIF

  706 CONTINUE

      WRITE (*,*) ' ENTER NAME OF INPUT FILE >'
      READ (*,701,ERR=706) INFILE
      DOTSUF='.'//SUFINP
      ISUF=INDEX(INFILE,DOTSUF)
      IF( ISUF.EQ.0 ) THEN
         JOBNAM = INFILE
      ELSE
         JOBNAM = INFILE(1:ISUF-1)
         REST=INFILE(ISUF+4:LEN(INFILE))
         ISUF=INDEX(REST,DOTSUF)
         IF( ISUF.NE.0 ) THEN
            WRITE(*,*) ' ERROR - INPUT FILE NAME HAS INVALID FORMAT'
            STOP
         ENDIF
      ENDIF

 707  CONTINUE

      WRITE (*,*) ' ENTER NAME OF OUTPUT FILE >'
      READ (*,701,ERR=707) OUTFIL
      CALL FILSUF(JOBNAM,FILBUF,0,SUFRES,IFAIL,DEBUG)
      IF(IFAIL.NE.0) THEN
         WRITE (*,*) ' ERROR CREATING OUTPUT FILE NAME'
         STOP
      ENDIF

      IFAIL=0
      CALL NEC2D( JOBNAM , FILBUF , IFAIL , DEBUG , MULTIF , LD , 
     &            IRESRV , NSMAX , LOADMX , NETMX , NORMF , MXCOUP , 
     &            MXCOU2 , JMAX , NSMAXX , NPMAX , IP , ICON1 , ICON2 , 
     &            ITAG , ICONX , IVQD , ISANT , IQDS , NCTAG , NCSEG , 
     &            ISEG1 , ISEG2 , NTYP , JCO , ISCON , IPCON , LDTYP , 
     &            LDTAG , LDTAGF , LDTAGT , IX , SALP , AIR , AII , 
     &            BIR , BII , CIR , CII , GAIN , X , Y  , Z , SI , BI , 
     &            ALP , BET , X2 , Y2 , Z2 ,  CAB , SAB , T1X , T1Y , 
     &            T1Z , T2X , T2Y , T2Z , XS , YS , ZS , Z2S , X11R , 
     &            X11I , X12R , X12I , X22R , X22I , AX , BX , CX , 
     &            ZLR , ZLI , ZLC , FNORM , XTEMP , YTEMP , ZTEMP , 
     &            SITEMP , BITEMP , CM , CUR , ZARRAY , D , VQD , 
     &            VSANT , VQDS , Y11A , Y12A , IPNT , NTEQA , NTSCA , 
     &            CMN , RHNT , RHS , VSRC , RHNX )
      IF(IFAIL.NE.0) THEN
         WRITE(*,*) ' CALL TO NEC2D FAILED WITH CODE ', IFAIL
      ENDIF

      STOP

  701 FORMAT (A)

      END
