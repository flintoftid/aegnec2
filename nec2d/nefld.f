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

      SUBROUTINE NEFLD( XOB , YOB , ZOB , EX , EY , EZ , LD , SALP ,
     &                  AIR , AII , BIR , BII , CIR , CII , CUR , N , 
     &                  M , ICON1 , ICON2 , X , Y , Z , SI , BI , 
     &                  CAB , SAB , T1X , T1Y , T1Z , T2X , T2Y , T2Z ,
     &                  IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* NEFLD computes the near field at specified points in space after   **
C* the structure currents have been computed.                         **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XOB                                                       **
C* INPUT  - YOB                                                       **
C* INPUT  - ZOB                                                       **
C* OUTPUT - EX                                                        **
C* OUTPUT - EY                                                        **
C* OUTPUT - EZ                                                        **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - AIR                                                       **
C* INPUT  - AII                                                       **
C* INPUT  - BIR                                                       **
C* INPUT  - BII                                                       **
C* INPUT  - CIR                                                       **
C* INPUT  - CII                                                       **
C* INPUT  - CUR                                                       **
C* INPUT  - N                                                         **
C* INPUT  - M                                                         **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - SI                                                        **
C* INPUT  - BI                                                        **
C* INPUT  - CAB                                                       **
C* INPUT  - SAB                                                       **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    B       CABJ    IND1    IND2    IPGND   S       SABJ   **
C*             SALPJ   T1XJ    T1YJ    T1ZJ    T2XJ    T2YJ    T2ZJ   **
C*             XJ      YJ      ZJ                                     **
C* uses value  B       CABJ    EXC     EXK     EXS     EYC     EYK    **
C*             EYS     EZC     EZK     EZS     IEXK    KSYMP   SABJ   **
C*             SALPJ   T1XJ    T1YJ    T1ZJ    T2XJ    T2YJ    T2ZJ   **
C*             XJ      YJ      ZJ                                     **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       EFLD    UNERE                                          **
C* called by   NFPAT   NHFLD                                          **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL EFLD , UNERE

C     Dummy arguments.
      INTEGER LD , N , M , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD)
      REAL*8 XOB, YOB ,ZOB 
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       BI(LD) , CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , 
     &       T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 EX , EY , EZ
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      INTEGER I , IIP , IPR , JC , JL
      REAL*8 XI , ZP , AAX
      COMPLEX*16 ACX , BCX , CCX

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
      INCLUDE 'gnd.inc'


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING NEFLD'
      ENDIF

      EX = (0.0D0,0.0D0)
      EY = (0.0D0,0.0D0)
      EZ = (0.0D0,0.0D0)
      AAX = 0.0D0
      IF ( N.EQ.0 ) GOTO 20
      DO 1 I = 1 , N
         XJ = XOB - X(I)
         YJ = YOB - Y(I)
         ZJ = ZOB - Z(I)
         ZP = CAB(I)*XJ + SAB(I)*YJ + SALP(I)*ZJ
         IF ( DABS(ZP).GT.0.5001D0*SI(I) ) GOTO 1
         ZP = XJ*XJ + YJ*YJ + ZJ*ZJ - ZP*ZP
         XJ = BI(I)
         IF ( ZP.GT.0.9D0*XJ*XJ ) GOTO 1
         AAX = XJ
         GOTO 2
    1 CONTINUE
    2 DO 19 I = 1 , N
         S = SI(I)
         B = BI(I)
         XJ = X(I)
         YJ = Y(I)
         ZJ = Z(I)
         CABJ = CAB(I)
         SABJ = SAB(I)
         SALPJ = SALP(I)
         IF ( IEXK.EQ.0 ) GOTO 18
         IPR = ICON1(I)
         IF ( IPR ) 3 , 8 , 4
    3    IPR = -IPR
         IF ( -ICON1(IPR).NE.I ) GOTO 9
         GOTO 6
    4    IF ( IPR.NE.I ) GOTO 5
         IF ( CABJ*CABJ+SABJ*SABJ.GT.1.D-8 ) GOTO 9
         GOTO 7
    5    IF ( ICON2(IPR).NE.I ) GOTO 9
    6    XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
         IF ( XI.LT.0.999999D0 ) GOTO 9
         IF ( DABS(BI(IPR)/B-1.0D0).GT.1.D-6 ) GOTO 9
    7    IND1 = 0
         GOTO 10
    8    IND1 = 1
         GOTO 10
    9    IND1 = 2
   10    IPR = ICON2(I)
         IF ( IPR ) 11 , 16 , 12
   11    IPR = -IPR
         IF ( -ICON2(IPR).NE.I ) GOTO 17
         GOTO 14
   12    IF ( IPR.NE.I ) GOTO 13
         IF ( CABJ*CABJ+SABJ*SABJ.GT.1.D-8 ) GOTO 17
         GOTO 15
   13    IF ( ICON1(IPR).NE.I ) GOTO 17
   14    XI = DABS(CABJ*CAB(IPR)+SABJ*SAB(IPR)+SALPJ*SALP(IPR))
         IF ( XI.LT.0.999999D0 ) GOTO 17
         IF ( DABS(BI(IPR)/B-1.0D0).GT.1.D-6 ) GOTO 17
   15    IND2 = 0
         GOTO 18
   16    IND2 = 1
         GOTO 18
   17    IND2 = 2
   18    CONTINUE
         CALL EFLD(XOB,YOB,ZOB,AAX,1,IFAIL)
         IF(IFAIL.NE.0) RETURN
         ACX = DCMPLX(AIR(I),AII(I))
         BCX = DCMPLX(BIR(I),BII(I))
         CCX = DCMPLX(CIR(I),CII(I))
         EX = EX + EXK*ACX + EXS*BCX + EXC*CCX
         EY = EY + EYK*ACX + EYS*BCX + EYC*CCX
         EZ = EZ + EZK*ACX + EZS*BCX + EZC*CCX
   19 CONTINUE
      IF ( M.EQ.0 ) RETURN
   20 JC = N
      JL = LD + 1
      DO 22 I = 1 , M
         JL = JL - 1
         S = BI(JL)
         XJ = X(JL)
         YJ = Y(JL)
         ZJ = Z(JL)
         T1XJ = T1X(JL)
         T1YJ = T1Y(JL)
         T1ZJ = T1Z(JL)
         T2XJ = T2X(JL)
         T2YJ = T2Y(JL)
         T2ZJ = T2Z(JL)
         JC = JC + 3
         ACX = T1XJ*CUR(JC-2) + T1YJ*CUR(JC-1) + T1ZJ*CUR(JC)
         BCX = T2XJ*CUR(JC-2) + T2YJ*CUR(JC-1) + T2ZJ*CUR(JC)
         DO 21 IIP = 1 , KSYMP
            IPGND = IIP
            CALL UNERE(XOB,YOB,ZOB)
            EX = EX + ACX*EXK + BCX*EXS
            EY = EY + ACX*EYK + BCX*EYS
            EZ = EZ + ACX*EZK + BCX*EZS
   21    CONTINUE
   22 CONTINUE
 
      RETURN
 
      END
