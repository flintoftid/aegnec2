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

      SUBROUTINE NHFLD( XOB , YOB , ZOB , HX , HY , HZ , LD , SALP ,
     &                  AIR , AII , BIR , BII , CIR , CII , CUR , N , 
     &                  M , ICON1 , ICON2 , X , Y , Z , SI , BI , CAB , 
     &                  SAB , T1X , T1Y , T1Z , T2X , T2Y , T2Z , 
     &                  IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* NHFLD computes the near H field at specified points in space after **
C* the structure currents have been computed.                         **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XOB                                                       **
C* INPUT  - YOB                                                       **
C* INPUT  - ZOB                                                       **
C* OUTPUT - HX                                                        **
C* OUTPUT - HY                                                        **
C* OUTPUT - HZ                                                        **
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
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
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
C* modifies    B       CABJ    S       SABJ    SALPJ   T1XJ    T1YJ   **
C*             T1ZJ    T2XJ    T2YJ    T2ZJ    XJ      YJ      ZJ     **
C* uses value  EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS     IPERF   T1XJ    T1YJ    T1ZJ    T2XJ   **
C*             T2YJ    T2ZJ    XJ      YJ      ZJ                     **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       HINTG   HSFLD   NEFLD                                  **
C* called by   NFPAT                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL HINTG , HSFLD , NEFLD

C     Duumy arguments.
      INTEGER LD , N , M , IFAIL , DEBUG 
      INTEGER ICON1(2*LD) , ICON2(2*LD)
      REAL*8 XOB , YOB , ZOB 
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       BI(LD) , CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , 
     &       T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 HX , HY , HZ
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      INTEGER I , JC , JL
      REAL*8 AAX , DELT , ZP
      COMPLEX*16 ACX , BCX , CCX , CON , EXPX , EXMX , EXPY , EXMY , 
     &           EXPZ , EXMZ , EYPX , EYMX , EYPY , EYMY , EYPZ , 
     &           EYMZ , EZPX , EZMX , EZPY , EZMY , EZPZ , EZMZ

C     Common storage.
      INCLUDE 'gnd.inc'
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING NHFLD'
      ENDIF

      IF ( IPERF.EQ.2 ) GOTO 6
      HX = (0.0D0,0.0D0)
      HY = (0.0D0,0.0D0)
      HZ = (0.0D0,0.0D0)
      AAX = 0.0D0
      IF ( N.EQ.0 ) GOTO 4
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
    2 DO 3 I = 1 , N
         S = SI(I)
         B = BI(I)
         XJ = X(I)
         YJ = Y(I)
         ZJ = Z(I)
         CABJ = CAB(I)
         SABJ = SAB(I)
         SALPJ = SALP(I)
         CALL HSFLD(XOB,YOB,ZOB,AAX)
         ACX = DCMPLX(AIR(I),AII(I))
         BCX = DCMPLX(BIR(I),BII(I))
         CCX = DCMPLX(CIR(I),CII(I))
         HX = HX + EXK*ACX + EXS*BCX + EXC*CCX
         HY = HY + EYK*ACX + EYS*BCX + EYC*CCX
         HZ = HZ + EZK*ACX + EZS*BCX + EZC*CCX
    3 CONTINUE
      IF ( M.EQ.0 ) RETURN
    4 JC = N
      JL = LD + 1
      DO 5 I = 1 , M
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
         CALL HINTG(XOB,YOB,ZOB,DEBUG)
         JC = JC + 3
         ACX = T1XJ*CUR(JC-2) + T1YJ*CUR(JC-1) + T1ZJ*CUR(JC)
         BCX = T2XJ*CUR(JC-2) + T2YJ*CUR(JC-1) + T2ZJ*CUR(JC)
         HX = HX + ACX*EXK + BCX*EXS
         HY = HY + ACX*EYK + BCX*EYS
         HZ = HZ + ACX*EZK + BCX*EZS
    5 CONTINUE
      RETURN

C     Get H by finite difference of E for sommerfeld ground
C     CON=j/(2*pi*eta).
C     DELT is the increment for getting central differences.

    6 DELT = 1.0D-3
      CON = (0.0D0,4.2246D-4)
      CALL NEFLD(XOB+DELT,YOB,ZOB,EXPX,EYPX,EZPX,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL NEFLD(XOB-DELT,YOB,ZOB,EXMX,EYMX,EZMX,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL NEFLD(XOB,YOB+DELT,ZOB,EXPY,EYPY,EZPY,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL NEFLD(XOB,YOB-DELT,ZOB,EXMY,EYMY,EZMY,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL NEFLD(XOB,YOB,ZOB+DELT,EXPZ,EYPZ,EZPZ,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      CALL NEFLD(XOB,YOB,ZOB-DELT,EXMZ,EYMZ,EZMZ,LD,SALP,AIR,AII,
     &           BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &           SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      HX = CON*(EZPY-EZMY-EYPZ+EYMZ)/(2.0D0*DELT)
      HY = CON*(EXPZ-EXMZ-EZPX+EZMX)/(2.0D0*DELT)
      HZ = CON*(EYPX-EYMX-EXPY+EXMY)/(2.0D0*DELT)
 
      RETURN
 
      END
