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
 
      SUBROUTINE MOVE( ROX , ROY , ROZ , XXS , YYS , ZZS , ITS , NRPT , 
     &                 ITGI , LD , SALP , IPSYM , N2 , N , NP , M1 , 
     &                 M2 , M , MP , ITAG , X , Y , Z , BI , X2 , Y2 , 
     &                 Z2 , T1X , T1Y , T1Z , T2X , T2Y , T2Z , IFAIL ,
     &                 DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine MOVE moves the structure with respect to its            **
C* coordinate system or reproduces structure in new positions.        **
C* Structure is rotated about X,Y,Z axes by ROX,ROY,ROZ               **
C* respectively, then shifted by XXS,YYS,ZZS.                         **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ROX                                                       **
C* INPUT  - ROY                                                       **
C* INPUT  - ROZ                                                       **
C* INPUT  - XXS                                                       **
C* INPUT  - YYS                                                       **
C* INPUT  - ZZS                                                       **
C* PASSED - ITS                                                       **
C* INPUT  - NRPT                                                      **
C* INPUT  - ITGI                                                      **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
C* OUTPUT - IPSYM                                                     **
C* INPUT  - N2                                                        **
C* OUTPUT - N                                                         **
C* OUTPUT - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M2                                                        **
C* OUTPUT - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - ITAG                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - BI                                                        **
C* OUTPUT - X2                                                        **
C* OUTPUT - Y2                                                        **
C* OUTPUT - Z2                                                        **
C* OUTPUT - T1X                                                       **
C* OUTPUT - T1Y                                                       **
C* OUTPUT - T1Z                                                       **
C* OUTPUT - T2X                                                       **
C* OUTPUT - T2Y                                                       **
C* OUTPUT - T2Z                                                       **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ISEGNO                                                 **
C* called by   DATAGN                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      INTEGER ISEGNO
      EXTERNAL ISEGNO

C     Dummy arguments.
      INTEGER ITS , NRPT , ITGI , LD , IPSYM , N2 , N , NP , M1 , M2 , 
     &        M , MP , IFAIL , DEBUG 
      INTEGER ITAG(2*LD)
      REAL*8 ROX , ROY , ROZ , XXS , YYS , ZZS
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , BI(LD) , X2(LD) , 
     &       Y2(LD) , Z2(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD)

C     Local variables.
      INTEGER I , I1 , II , IR , IX , K , KR , LDI , NRP
      REAL*8 CPH , CPS , CTH , SPH , SPS , STH , XI , XX , XY , XZ , 
     &       YI , YX , YY , YZ , ZI , ZX , ZY , ZZ


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING MOVE'
      ENDIF

      IF ( DABS(ROX)+DABS(ROY).GT.1.D-10 ) IPSYM = IPSYM*3
      SPS = DSIN(ROX)
      CPS = DCOS(ROX)
      STH = DSIN(ROY)
      CTH = DCOS(ROY)
      SPH = DSIN(ROZ)
      CPH = DCOS(ROZ)
      XX = CPH*CTH
      XY = CPH*STH*SPS - SPH*CPS
      XZ = CPH*STH*CPS + SPH*SPS
      YX = SPH*CTH
      YY = SPH*STH*SPS + CPH*CPS
      YZ = SPH*STH*CPS - CPH*SPS
      ZX = -STH
      ZY = CTH*SPS
      ZZ = CTH*CPS
      NRP = NRPT
      IF ( NRPT.EQ.0 ) NRP = 1
      IX = 1
      IF ( N.LT.N2 ) GOTO 3
      I1 = ISEGNO(ITS,1,LD,ITAG,N,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      IF ( I1.LT.N2 ) I1 = N2
      IX = I1
C     IDF. Fix array overflow bug.
      IF((M+N+NRPT*(N-IX+1)).GT.LD) THEN
         IFAIL = 47
         RETURN
      ENDIF
C     End IDF.
      K = N
      IF ( NRPT.EQ.0 ) K = I1 - 1
      DO 2 IR = 1 , NRP
         DO 1 I = I1 , N
            K = K + 1
            XI = X(I)
            YI = Y(I)
            ZI = Z(I)
            X(K) = XI*XX + YI*XY + ZI*XZ + XXS
            Y(K) = XI*YX + YI*YY + ZI*YZ + YYS
            Z(K) = XI*ZX + YI*ZY + ZI*ZZ + ZZS
            XI = X2(I)
            YI = Y2(I)
            ZI = Z2(I)
            X2(K) = XI*XX + YI*XY + ZI*XZ + XXS
            Y2(K) = XI*YX + YI*YY + ZI*YZ + YYS
            Z2(K) = XI*ZX + YI*ZY + ZI*ZZ + ZZS
            BI(K) = BI(I)
            ITAG(K) = ITAG(I)
            IF ( ITAG(I).NE.0 ) ITAG(K) = ITAG(I) + ITGI
    1    CONTINUE
         I1 = N + 1
         N = K
    2 CONTINUE
    3 IF ( M.LT.M2 ) GOTO 6
      I1 = M2
C     IDF. Fix array overflow bug.
      IF((M+N+NRPT*(M-I1+1)).GT.LD) THEN
         IFAIL = 47
         RETURN
      ENDIF
C     End IDF.
      K = M
      LDI = LD + 1
      IF ( NRPT.EQ.0 ) K = M1
      DO 5 II = 1 , NRP
         DO 4 I = I1 , M
            K = K + 1
            IR = LDI - I
            KR = LDI - K
            XI = X(IR)
            YI = Y(IR)
            ZI = Z(IR)
            X(KR) = XI*XX + YI*XY + ZI*XZ + XXS
            Y(KR) = XI*YX + YI*YY + ZI*YZ + YYS
            Z(KR) = XI*ZX + YI*ZY + ZI*ZZ + ZZS
            XI = T1X(IR)
            YI = T1Y(IR)
            ZI = T1Z(IR)
            T1X(KR) = XI*XX + YI*XY + ZI*XZ
            T1Y(KR) = XI*YX + YI*YY + ZI*YZ
            T1Z(KR) = XI*ZX + YI*ZY + ZI*ZZ
            XI = T2X(IR)
            YI = T2Y(IR)
            ZI = T2Z(IR)
            T2X(KR) = XI*XX + YI*XY + ZI*XZ
            T2Y(KR) = XI*YX + YI*YY + ZI*YZ
            T2Z(KR) = XI*ZX + YI*ZY + ZI*ZZ
            SALP(KR) = SALP(IR)
            BI(KR) = BI(IR)
    4    CONTINUE
         I1 = M + 1
         M = K
    5 CONTINUE
    6 IF ( (NRPT.EQ.0) .AND. (IX.EQ.1) ) RETURN
      NP = N
      MP = M
      IPSYM = 0
 
      RETURN
 
      END
