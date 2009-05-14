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

#include <config.h>

C*--------------------------------------------------------------------**

      SUBROUTINE NETWK( CCM , CMB , CMC , CMD , IIP , EINC , FFMHZ , 
     &                  LD , AIR , AII , BIR , BII , CIR , CII , 
     &                  D , N1 , N , NP , M1 , M , MP , ICON1 , 
     &                  ICON2 , ITAG , WLAM , SI , BI , T1X , T1Y , 
     &                  T1Z , T2X , T2Y , T2Z , NSMAX , NVQD , NSANT , 
     &                  NQDS , IVQD , ISANT , IQDS , VQD , VSANT , 
     &                  VQDS , NETMX , NEQ , NEQ2 , NONET , NTSOL , 
     &                  NPRINT , MASYM , ISEG1 , ISEG2 , NTYP , PIN , 
     &                  PNLS , X11R , X11I , X12R , X12I , X22R , X22I , 
     &                  ZPED , JMAX , NSMAXX , JSNO , NSCON , NPCON , 
     &                  JCO , ISCON , AX , BX , CX , MULTIF , IPNT , 
     &                  NTEQA , NTSCA , CMN , RHNT , RHS , VSRC , RHNX,
     &                  IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine NETWK solves for the voltages and currents at the posts **
C* of non-radiating networks that are part of the antenna. NETWK is   **
C* also involved in the solution when no non-radiating networks are   **
C* present and computes the relative driving point matrix asymmetry   **
C* when this option is requested.                                     **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - CCM                                                       **
C* PASSED - CMB                                                       **
C* PASSED - CMC                                                       **
C* PASSED - CMD                                                       **
C* PASSED - IIP                                                       **
C* OUTPUT - EINC                                                      **
C* INPUT  - FFMHZ                                                     **
C* INPUT  - LD                                                        **
C* INPUT  - AIR                                                       **
C* INPUT  - AII                                                       **
C* INPUT  - BIR                                                       **
C* INPUT  - BII                                                       **
C* INPUT  - CIR                                                       **
C* INPUT  - CII                                                       **
C* PASSED - D                                                         **
C* PASSED - N1                                                        **
C* PASSED - N                                                         **
C* PASSED - NP                                                        **
C* PASSED - M1                                                        **
C* PASSED - M                                                         **
C* PASSED - MP                                                        **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - ITAG                                                      **
C* INPUT  - WLAM                                                      **
C* INPUT  - SI                                                        **
C* PASSED - BI                                                        **
C* PASSED - T1X                                                       **
C* PASSED - T1Y                                                       **
C* PASSED - T1Z                                                       **
C* PASSED - T2X                                                       **
C* PASSED - T2Y                                                       **
C* PASSED - T2Z                                                       **
C* INPUT  - NSMAX                                                     **
C* INPUT  - NVQD                                                      **
C* INPUT  - NSANT                                                     **
C* PASSED - NQDS                                                      **
C* INPUT  - IVQD                                                      **
C* INPUT  - ISANT                                                     **
C* PASSED - IQDS                                                      **
C* INPUT  - VQD                                                       **
C* INPUT  - VSANT                                                     **
C* PASSED - VQDS                                                      **
C* INPUT  - NETMX                                                     **
C* INPUT  - NEQ                                                       **
C* INPUT  - NEQ2                                                      **
C* INPUT  - NONET                                                     **
C* INPUT  - NTSOL                                                     **
C* INPUT  - NPRINT                                                    **
C* INPUT  - MASYM                                                     **
C* INPUT  - ISEG1                                                     **
C* INPUT  - ISEG2                                                     **
C* INPUT  - NTYP                                                      **
C* OUTPUT - PIN                                                       **
C* OUTPUT - PNLS                                                      **
C* INPUT  - X11R                                                      **
C* INPUT  - X11I                                                      **
C* INPUT  - X12R                                                      **
C* INPUT  - X12I                                                      **
C* INPUT  - X22R                                                      **
C* INPUT  - X22I                                                      **
C* OUTPUT - ZPED                                                      **
C* INPUT  - JMAX                                                      **
C* INPUT  - NSMAXX                                                    **
C* PASSED - JSNO                                                      **
C* PASSED - NSCON                                                     **
C* PASSED - NPCON                                                     **
C* PASSED - JCO                                                       **
C* PASSED - ISCON                                                     **
C* PASSED - AX                                                        **
C* PASSED - BX                                                        **
C* PASSED - CX                                                        **
C* INPUT  - MULTIF                                                    **
C* OUTPUT - IPNT                                                      **
C* OUTPUT - NTEQA                                                     **
C* OUTPUT - NTSCA                                                     **
C* OUTPUT - CMN                                                       **
C* OUTPUT - RHNT                                                      **
C* OUTPUT - RHS                                                       **
C* OUTPUT - VSRC                                                      **
C* OUTPUT - RHNX                                                      **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       CABC    FACTR   SOLGF   SOLVE                          **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL  CABC , FACTR , SOLGF , SOLVE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER LD , N1 , N , NP , M1 , M , MP , NSMAX , NVQD , 
     &        NSANT , NQDS , NETMX , NEQ , NEQ2 , NONET , NTSOL , 
     &        NPRINT , MASYM , JMAX , NSMAXX , JSNO , NSCON , NPCON ,
     &        MULTIF , IFAIL , DEBUG
      INTEGER IIP(1) , ICON1(2*LD) , ICON2(2*LD) , ITAG(2*LD) , 
     &        IVQD(NSMAX) , ISANT(NSMAX) , IQDS(NSMAX) , ISEG1(NETMX) , 
     &        ISEG2(NETMX) , NTYP(NETMX) , JCO(JMAX) , ISCON(NSMAXX) , 
     &        IPNT(NETMX) , NTEQA(NETMX) , NTSCA(NETMX)
      REAL*8 FFMHZ , WLAM , PIN , PNLS 
      REAL*8 AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , SI(LD) , BI(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , X11R(NETMX) , X11I(NETMX) , 
     &       X12R(NETMX) , X12I(NETMX) , X22R(NETMX) , X22I(NETMX) , 
     &       AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 ZPED 
      COMPLEX*16 CCM(1) , CMB(1) , CMC(1) , CMD(1) , EINC(1) , 
     &           D(2*LD) , VQD(NSMAX) , VSANT(NSMAX) , VQDS(NSMAX) ,
     &           CMN(NETMX,NETMX) , RHNT(NETMX) , RHS(3*LD) , 
     &           VSRC(NETMX) , RHNX(NETMX)

C     Local variables.
      INTEGER I , IROW1 , IROW2 , ISC1 , ISC2 , J , NDIMN , NDIMNP , 
     &        NEQT , NEQZ2 , NSEG1 , NSEG2 , NTEQ , NTSC 
      REAL*8 ASA , ASM , PWR , TP , Y11I , Y11R , Y12I , Y12R , 
     &       Y22I , Y22R 
      COMPLEX*16 YMIT , VLT , CUX

C     Data initialisation.
      DATA TP /6.283185307179586476925286D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING NETWK'
      ENDIF
 
      NDIMN = NETMX
      NDIMNP = NDIMN + 1

      NEQZ2 = NEQ2
      IF ( NEQZ2.EQ.0 ) NEQZ2 = 1
      PIN = 0.0D0
      PNLS = 0.0D0
      NEQT = NEQ + NEQ2
      IF ( NTSOL.NE.0 ) GOTO 46
      IF ( MASYM.EQ.0 ) GOTO 17

C     Compute relative matrix asymmetry.

      IROW1 = 0
      IF ( NONET.EQ.0 ) GOTO 6
      DO 5 I = 1 , NONET
         NSEG1 = ISEG1(I)
         DO 4 ISC1 = 1 , 2
            IF ( IROW1.EQ.0 ) GOTO 2
            DO 1 J = 1 , IROW1
               IF ( NSEG1.EQ.IPNT(J) ) GOTO 3
    1       CONTINUE
    2       IROW1 = IROW1 + 1
            IPNT(IROW1) = NSEG1
    3       NSEG1 = ISEG2(I)
    4    CONTINUE
    5 CONTINUE
    6 IF ( NSANT.EQ.0 ) GOTO 10
      DO 9 I = 1 , NSANT
         NSEG1 = ISANT(I)
         IF ( IROW1.EQ.0 ) GOTO 8
         DO 7 J = 1 , IROW1
            IF ( NSEG1.EQ.IPNT(J) ) GOTO 9
    7    CONTINUE
    8    IROW1 = IROW1 + 1
         IPNT(IROW1) = NSEG1
    9 CONTINUE
   10 IF ( IROW1.LT.NDIMNP ) GOTO 11
      WRITE (CHRSLT,63)
      IFAIL=27
      RETURN
   11 IF ( IROW1.LT.2 ) GOTO 17
      DO 14 I = 1 , IROW1
         ISC1 = IPNT(I)
         ASM = SI(ISC1)
         DO 12 J = 1 , NEQT
            RHS(J) = DCMPLX(0.0D0,0.0D0)
   12    CONTINUE
         RHS(ISC1) = DCMPLX(1.0D0,0.0D0)
         CALL SOLGF(CCM,CMB,CMC,CMD,RHS,IIP,NP,N1,N,MP,M1,M,NEQ,NEQ2,
     &              NEQZ2,LD,D,NSMAXX,NSCON,NPCON,ISCON,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         CALL CABC(RHS,LD,AIR,AII,BIR,BII,CIR,CII,N,M,ICON1,ICON2,WLAM,
     &             SI,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NQDS,IQDS,VQDS,
     &             JMAX,JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         DO 13 J = 1 , IROW1
            ISC1 = IPNT(J)
            CMN(J,I) = RHS(ISC1)/ASM
   13    CONTINUE
   14 CONTINUE
      ASM = 0.0D0
      ASA = 0.0D0
      DO 16 I = 2 , IROW1
         ISC1 = I - 1
         DO 15 J = 1 , ISC1
            CUX = CMN(I,J)
            PWR = CDABS((CUX-CMN(J,I))/CUX)
            ASA = ASA + PWR*PWR
            IF ( PWR.LT.ASM ) GOTO 15
            ASM = PWR
            NTEQ = IPNT(I)
            NTSC = IPNT(J)
   15    CONTINUE
   16 CONTINUE
      ASA = DSQRT(ASA*2.0D0/DBLE(IROW1*(IROW1-1)))
      WRITE (CHRSLT,62) ASM , NTEQ , NTSC , ASA
   17 IF ( NONET.EQ.0 ) GOTO 52

C     Solution of network equations.

      DO 19 I = 1 , NDIMN
         RHNX(I) = DCMPLX(0.0D0,0.0D0)
         DO 18 J = 1 , NDIMN
            CMN(I,J) = DCMPLX(0.0D0,0.0D0)
   18    CONTINUE
   19 CONTINUE
      NTEQ = 0
      NTSC = 0

C     Sort network and source data and assign equation numbers to
C     segments.

      DO 42 J = 1 , NONET
         NSEG1 = ISEG1(J)
         NSEG2 = ISEG2(J)
         IF ( NTYP(J).GT.1 ) GOTO 20
         Y11R = X11R(J)
         Y11I = X11I(J)
         Y12R = X12R(J)
         Y12I = X12I(J)
         Y22R = X22R(J)
         Y22I = X22I(J)
         GOTO 21
   20    Y22R = TP*X11I(J)/WLAM
         Y12R = 0.0D0
         Y12I = 1.0D0/(X11R(J)*DSIN(Y22R))
         Y11R = X12R(J)
         Y11I = -Y12I*DCOS(Y22R)
         Y22R = X22R(J)
         Y22I = Y11I + X22I(J)
         Y11I = Y11I + X12I(J)
         IF ( NTYP(J).EQ.2 ) GOTO 21
         Y12R = -Y12R
         Y12I = -Y12I
   21    IF ( NSANT.EQ.0 ) GOTO 23
         DO 22 I = 1 , NSANT
            IF ( NSEG1.NE.ISANT(I) ) GOTO 22
            ISC1 = I
            GOTO 26
   22    CONTINUE
   23    ISC1 = 0
         IF ( NTEQ.EQ.0 ) GOTO 25
         DO 24 I = 1 , NTEQ
            IF ( NSEG1.NE.NTEQA(I) ) GOTO 24
            IROW1 = I
            GOTO 29
   24    CONTINUE
   25    NTEQ = NTEQ + 1
         IROW1 = NTEQ
         NTEQA(NTEQ) = NSEG1
         GOTO 29
   26    IF ( NTSC.EQ.0 ) GOTO 28
         DO 27 I = 1 , NTSC
            IF ( NSEG1.NE.NTSCA(I) ) GOTO 27
            IROW1 = NDIMNP - I
            GOTO 29
   27    CONTINUE
   28    NTSC = NTSC + 1
         IROW1 = NDIMNP - NTSC
         NTSCA(NTSC) = NSEG1
         VSRC(NTSC) = VSANT(ISC1)
   29    IF ( NSANT.EQ.0 ) GOTO 31
         DO 30 I = 1 , NSANT
            IF ( NSEG2.NE.ISANT(I) ) GOTO 30
            ISC2 = I
            GOTO 34
   30    CONTINUE
   31    ISC2 = 0
         IF ( NTEQ.EQ.0 ) GOTO 33
         DO 32 I = 1 , NTEQ
            IF ( NSEG2.NE.NTEQA(I) ) GOTO 32
            IROW2 = I
            GOTO 37
   32    CONTINUE
   33    NTEQ = NTEQ + 1
         IROW2 = NTEQ
         NTEQA(NTEQ) = NSEG2
         GOTO 37
   34    IF ( NTSC.EQ.0 ) GOTO 36
         DO 35 I = 1 , NTSC
            IF ( NSEG2.NE.NTSCA(I) ) GOTO 35
            IROW2 = NDIMNP - I
            GOTO 37
   35    CONTINUE
   36    NTSC = NTSC + 1
         IROW2 = NDIMNP - NTSC
         NTSCA(NTSC) = NSEG2
         VSRC(NTSC) = VSANT(ISC2)
   37    IF ( NTSC+NTEQ.LT.NDIMNP ) GOTO 38
         WRITE (CHRSLT,63)
         IFAIL=27
         RETURN

C     Fill network equation matrix and right hand side vector with
C     network short-circuit admittance matrix coefficients.

   38    IF ( ISC1.NE.0 ) GOTO 39
         CMN(IROW1,IROW1) = CMN(IROW1,IROW1) - DCMPLX(Y11R,Y11I)
     &                      *SI(NSEG1)
         CMN(IROW1,IROW2) = CMN(IROW1,IROW2) - DCMPLX(Y12R,Y12I)
     &                      *SI(NSEG1)
         GOTO 40
   39    RHNX(IROW1) = RHNX(IROW1) + DCMPLX(Y11R,Y11I)*VSANT(ISC1)/WLAM
         RHNX(IROW2) = RHNX(IROW2) + DCMPLX(Y12R,Y12I)*VSANT(ISC1)/WLAM
   40    IF ( ISC2.NE.0 ) GOTO 41
         CMN(IROW2,IROW2) = CMN(IROW2,IROW2) - DCMPLX(Y22R,Y22I)
     &                      *SI(NSEG2)
         CMN(IROW2,IROW1) = CMN(IROW2,IROW1) - DCMPLX(Y12R,Y12I)
     &                      *SI(NSEG2)
         GOTO 42
   41    RHNX(IROW1) = RHNX(IROW1) + DCMPLX(Y12R,Y12I)*VSANT(ISC2)/WLAM
         RHNX(IROW2) = RHNX(IROW2) + DCMPLX(Y22R,Y22I)*VSANT(ISC2)/WLAM
   42 CONTINUE

C     Add interaction matrix admittance elements to network equation
C     matrix.

      DO 45 I = 1 , NTEQ
         DO 43 J = 1 , NEQT
            RHS(J) = DCMPLX(0.0D0,0.0D0)
   43    CONTINUE
         IROW1 = NTEQA(I)
         RHS(IROW1) = DCMPLX(1.0D0,0.0D0)
         CALL SOLGF(CCM,CMB,CMC,CMD,RHS,IIP,NP,N1,N,MP,M1,M,NEQ,NEQ2,
     &              NEQZ2,LD,D,NSMAXX,NSCON,NPCON,ISCON,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         CALL CABC(RHS,LD,AIR,AII,BIR,BII,CIR,CII,N,M,ICON1,ICON2,WLAM,
     &             SI,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NQDS,IQDS,VQDS,
     &             JMAX,JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         DO 44 J = 1 , NTEQ
            IROW1 = NTEQA(J)
            CMN(I,J) = CMN(I,J) + RHS(IROW1)
   44    CONTINUE
   45 CONTINUE

C     Factor network equation matrix.

      CALL FACTR(NTEQ,CMN,IPNT,NDIMN,LD,D,DEBUG)

C     Add to network equation right hand side the terms due to element
C     interactions.

   46 IF ( NONET.EQ.0 ) GOTO 52
      DO 47 I = 1 , NEQT
         RHS(I) = EINC(I)
   47 CONTINUE
      CALL SOLGF(CCM,CMB,CMC,CMD,RHS,IIP,NP,N1,N,MP,M1,M,NEQ,NEQ2,NEQZ2,
     &           LD,D,NSMAXX,NSCON,NPCON,ISCON,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL CABC(RHS,LD,AIR,AII,BIR,BII,CIR,CII,N,M,ICON1,ICON2,WLAM,SI,
     &          BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NQDS,IQDS,VQDS,JMAX,
     &          JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      DO 48 I = 1 , NTEQ
         IROW1 = NTEQA(I)
         RHNT(I) = RHNX(I) + RHS(IROW1)
   48 CONTINUE

C     Solve network equations.

      CALL SOLVE(NTEQ,CMN,IPNT,RHNT,NDIMN,LD,D,DEBUG)

C     Add fields due to network voltages to electric fields applied to
C     structure and solve for induced current.

      DO 49 I = 1 , NTEQ
         IROW1 = NTEQA(I)
         EINC(IROW1) = EINC(IROW1) - RHNT(I)
   49 CONTINUE
      CALL SOLGF(CCM,CMB,CMC,CMD,EINC,IIP,NP,N1,N,MP,M1,M,NEQ,NEQ2,
     &           NEQZ2,LD,D,NSMAXX,NSCON,NPCON,ISCON,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL CABC(EINC,LD,AIR,AII,BIR,BII,CIR,CII,N,M,ICON1,ICON2,WLAM,SI,
     &          BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NQDS,IQDS,VQDS,JMAX,
     &          JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IF ( NPRINT.EQ.0 ) WRITE (CHRSLT,65)
      IF ( NPRINT.EQ.0 ) WRITE (CHRSLT,64)
      DO 50 I = 1 , NTEQ
         IROW1 = NTEQA(I)
         VLT = RHNT(I)*SI(IROW1)*WLAM
         CUX = EINC(IROW1)*WLAM
         YMIT = CUX/VLT
         ZPED = VLT/CUX
         IROW2 = ITAG(IROW1)
         PWR = 0.5D0*DBLE(VLT*DCONJG(CUX))
         PNLS = PNLS - PWR
         IF ( NPRINT.EQ.0 ) THEN
            WRITE (CHRSLT,66) IROW2 , IROW1 , VLT , 
     &                                        CUX , ZPED , YMIT , PWR
            IF ( MULTIF.GT.0 ) THEN
               WRITE(CHIMPD,641) FFMHZ, ZPED
#ifdef HAVE_FLUSH
	       CALL FLUSH(CHIMPD)
#endif
            ENDIF
         ENDIF
   50 CONTINUE
      IF ( NTSC.EQ.0 ) GOTO 53
      DO 51 I = 1 , NTSC
         IROW1 = NTSCA(I)
         VLT = VSRC(I)
         CUX = EINC(IROW1)*WLAM
         YMIT = CUX/VLT
         ZPED = VLT/CUX
         IROW2 = ITAG(IROW1)
         PWR = 0.5D0*DBLE(VLT*DCONJG(CUX))
         PNLS = PNLS - PWR
         IF ( NPRINT.EQ.0 ) THEN
            WRITE (CHRSLT,66) IROW2 , IROW1 , VLT , 
     &                                        CUX , ZPED , YMIT , PWR
            IF ( MULTIF.GT.0 ) THEN
               WRITE(CHIMPD,641) FFMHZ, ZPED
#ifdef HAVE_FLUSH
               CALL FLUSH(CHIMPD)
#endif
            ENDIF
         ENDIF
   51 CONTINUE
      GOTO 53

C     Solve for currents when no networks are present.

   52 CALL SOLGF(CCM,CMB,CMC,CMD,EINC,IIP,NP,N1,N,MP,M1,M,NEQ,NEQ2,
     &           NEQZ2,LD,D,NSMAXX,NSCON,NPCON,ISCON,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL CABC(EINC,LD,AIR,AII,BIR,BII,CIR,CII,N,M,ICON1,ICON2,WLAM,SI,
     &          BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,NSMAX,NQDS,IQDS,VQDS,JMAX,
     &          JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      NTSC = 0
   53 IF ( NSANT+NVQD.EQ.0 ) RETURN
      WRITE (CHRSLT,67)
      WRITE (CHRSLT,64)
      IF ( NSANT.EQ.0 ) GOTO 60
      DO 59 I = 1 , NSANT
         ISC1 = ISANT(I)
         VLT = VSANT(I)
         IF ( NTSC.EQ.0 ) GOTO 55
         DO 54 J = 1 , NTSC
            IF ( NTSCA(J).EQ.ISC1 ) GOTO 56
   54    CONTINUE
   55    CUX = EINC(ISC1)*WLAM
         IROW1 = 0
         GOTO 58
   56    IROW1 = NDIMNP - J
         CUX = RHNX(IROW1)
         DO 57 J = 1 , NTEQ
            CUX = CUX - CMN(J,IROW1)*RHNT(J)
   57    CONTINUE
         CUX = (EINC(ISC1)+CUX)*WLAM
   58    YMIT = CUX/VLT
         ZPED = VLT/CUX
         PWR = 0.5D0*DBLE(VLT*DCONJG(CUX))
         PIN = PIN + PWR
         IF ( IROW1.NE.0 ) PNLS = PNLS + PWR
         IROW2 = ITAG(ISC1)
         WRITE (CHRSLT,66) IROW2 , ISC1 , VLT , CUX , ZPED , YMIT , PWR
         IF ( MULTIF.GT.0 ) THEN
            WRITE(CHIMPD,641) FFMHZ, ZPED
#ifdef HAVE_FLUSH
            CALL FLUSH(CHIMPD)
#endif

         ENDIF
   59 CONTINUE
   60 IF ( NVQD.EQ.0 ) RETURN
      DO 61 I = 1 , NVQD
         ISC1 = IVQD(I)
         VLT = VQD(I)
         CUX = DCMPLX(AIR(ISC1),AII(ISC1))
         YMIT = DCMPLX(BIR(ISC1),BII(ISC1))
         ZPED = DCMPLX(CIR(ISC1),CII(ISC1))
         PWR = SI(ISC1)*TP*0.5D0
         CUX = (CUX-YMIT*DSIN(PWR)+ZPED*DCOS(PWR))*WLAM
         YMIT = CUX/VLT
         ZPED = VLT/CUX
         PWR = 0.5D0*DBLE(VLT*DCONJG(CUX))
         PIN = PIN + PWR
         IROW2 = ITAG(ISC1)
         WRITE (CHRSLT,68) IROW2 , ISC1 , VLT , CUX , ZPED , YMIT , PWR
         IF ( MULTIF.GT.0 ) THEN
            WRITE(CHIMPD,641) FFMHZ, ZPED
#ifdef HAVE_FLUSH
            CALL FLUSH(CHIMPD)
#endif
         ENDIF
   61 CONTINUE
 
      RETURN
 
   62 FORMAT (///,3X,'MAXIMUM RELATIVE ASYMMETRY OF THE DRIVING POINT',
     &        ' ADMITTANCE MATRIX IS',1P,E10.3,' FOR SEGMENTS',I5,
     &        ' AND',I5,/,3X,'RMS RELATIVE ASYMMETRY IS',E10.3)
   63 FORMAT (1X,'ERROR - - NETWORK ARRAY DIMENSIONS TOO SMALL')
   64 FORMAT (/,3X,'TAG',3X,'SEG.',4X,'VOLTAGE (VOLTS)',9X,
     &        'CURRENT (AMPS)',9X,'IMPEDANCE (OHMS)',8X,
     &        'ADMITTANCE (MHOS)',6X,'POWER',/,3X,'NO.',3X,'NO.',4X,
     &        'REAL',8X,'IMAG.',3(7X,'REAL',8X,'IMAG.'),5X,'(WATTS)')
   65 FORMAT (///,27X,
     &'- - - STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS - - 
     &-')
   66 FORMAT (2(1X,I5),1P,9E12.5)
   67 FORMAT (///,42X,'- - - ANTENNA INPUT PARAMETERS - - -')
   68 FORMAT (1X,I5,' *',I4,1P,9E12.5)
  641 FORMAT (2X,E12.5,4X,E12.5,4X,E12.5)
      END
