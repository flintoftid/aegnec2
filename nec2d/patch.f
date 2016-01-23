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

      SUBROUTINE PATCH( NX , NY , X1 , Y1 , Z1 , XXX2 , YYY2 , ZZZ2 , 
     &                  X3 , Y3 , Z3 , X4 , Y4 , Z4 , LD , SALP , 
     &                  IPSYM , N , NP , M , MP , X , Y, Z , BI , T1X , 
     &                  T1Y , T1Z , T2X , T2Y , T2Z , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* PATCH generates and modifies patch geometry data.                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NX                                                        **
C* INPUT  - NY                                                        **
C* INPUT  - X1                                                        **
C* INPUT  - Y1                                                        **
C* INPUT  - Z1                                                        **
C* INPUT  - XXX2                                                      **
C* INPUT  - YYY2                                                      **
C* INPUT  - ZZZ2                                                      **
C* INPUT  - X3                                                        **
C* INPUT  - Y3                                                        **
C* INPUT  - Z3                                                        **
C* INPUT  - X4                                                        **
C* INPUT  - Y4                                                        **
C* INPUT  - Z4                                                        **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
C* OUTPUT - IPSYM                                                     **
C* INPUT  - N                                                         **
C* OUTPUT - NP                                                        **
C* OUTPUT - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - BI                                                        **
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
C* calls       ** NOTHING **                                          **
C* called by   DATAGN                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NX , NY , LD , IPSYM , N , NP , M , MP , IFAIL , DEBUG
      REAL*8 X1 , Y1 , Z1 , XXX2 ,YYY2 , ZZZ2 , X3 ,Y3 , Z3 , X4 , 
     &       Y4 , Z4 
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , BI(LD) , T1X(LD) , 
     &       T1Y(LD) , T1Z(LD) , T2X(LD) , T2Y(LD) , T2Z(LD)

C     Local variables.
      INTEGER IX , IY , NTP , MI
      REAL*8 S1X , S1Y , S1Z , S2X , S2Y , S2Z , SALPN , XA , XN2 , 
     &       XNV , XXS , XST , XT , YN2 , YNV , YYS , YT , ZN2 , ZNV , 
     &       ZZS , ZT


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING PATCH'
      ENDIF

C     New patches.  For NX=0, NY=1,2,3,4 patch is (respectively)
C     arbitrary, rectagular, triangular, or quadrilateral.
C     For NX and NY .GT. 0 a rectangular surface is produced with
C     NX by NY rectangular patches.
      M = M + 1
      MI = LD + 1 - M
      NTP = NY
      IF ( NX.GT.0 ) NTP = 2
      IF ( NTP.GT.1 ) GOTO 2
      X(MI) = X1
      Y(MI) = Y1
      Z(MI) = Z1
      BI(MI) = ZZZ2
      ZNV = DCOS(XXX2)
      XNV = ZNV*DCOS(YYY2)
      YNV = ZNV*DSIN(YYY2)
      ZNV = DSIN(XXX2)
      XA = DSQRT(XNV*XNV+YNV*YNV)
      IF ( XA.LT.1.D-6 ) GOTO 1
      T1X(MI) = -YNV/XA
      T1Y(MI) = XNV/XA
      T1Z(MI) = 0.0D0
      GOTO 6
    1 T1X(MI) = 1.0D0
      T1Y(MI) = 0.0D0
      T1Z(MI) = 0.0D0
      GOTO 6
    2 S1X = XXX2 - X1
      S1Y = YYY2 - Y1
      S1Z = ZZZ2 - Z1
      S2X = X3 - XXX2
      S2Y = Y3 - YYY2
      S2Z = Z3 - ZZZ2
      IF ( NX.EQ.0 ) GOTO 3
      S1X = S1X/NX
      S1Y = S1Y/NX
      S1Z = S1Z/NX
      S2X = S2X/NY
      S2Y = S2Y/NY
      S2Z = S2Z/NY
    3 XNV = S1Y*S2Z - S1Z*S2Y
      YNV = S1Z*S2X - S1X*S2Z
      ZNV = S1X*S2Y - S1Y*S2X
      XA = DSQRT(XNV*XNV+YNV*YNV+ZNV*ZNV)
      XNV = XNV/XA
      YNV = YNV/XA
      ZNV = ZNV/XA
      XST = DSQRT(S1X*S1X+S1Y*S1Y+S1Z*S1Z)
      T1X(MI) = S1X/XST
      T1Y(MI) = S1Y/XST
      T1Z(MI) = S1Z/XST
      IF ( NTP.GT.2 ) GOTO 4
      X(MI) = X1 + 0.5D0*(S1X+S2X)
      Y(MI) = Y1 + 0.5D0*(S1Y+S2Y)
      Z(MI) = Z1 + 0.5D0*(S1Z+S2Z)
      BI(MI) = XA
      GOTO 6
    4 IF ( NTP.EQ.4 ) GOTO 5
      X(MI) = (X1+XXX2+X3)/3.0D0
      Y(MI) = (Y1+YYY2+Y3)/3.0D0
      Z(MI) = (Z1+ZZZ2+Z3)/3.0D0
      BI(MI) = 0.5D0*XA
      GOTO 6
    5 S1X = X3 - X1
      S1Y = Y3 - Y1
      S1Z = Z3 - Z1
      S2X = X4 - X1
      S2Y = Y4 - Y1
      S2Z = Z4 - Z1
      XN2 = S1Y*S2Z - S1Z*S2Y
      YN2 = S1Z*S2X - S1X*S2Z
      ZN2 = S1X*S2Y - S1Y*S2X
      XST = DSQRT(XN2*XN2+YN2*YN2+ZN2*ZN2)
      SALPN = 1.0D0/(3.0D0*(XA+XST))
      X(MI) = (XA*(X1+XXX2+X3)+XST*(X1+X3+X4))*SALPN
      Y(MI) = (XA*(Y1+YYY2+Y3)+XST*(Y1+Y3+Y4))*SALPN
      Z(MI) = (XA*(Z1+ZZZ2+Z3)+XST*(Z1+Z3+Z4))*SALPN
      BI(MI) = 0.5D0*(XA+XST)
      S1X = (XNV*XN2+YNV*YN2+ZNV*ZN2)/XST
      IF ( S1X.GT.0.9998D0 ) GOTO 6
      WRITE (CHRSLT,14)
      IFAIL=29
      RETURN
    6 T2X(MI) = YNV*T1Z(MI) - ZNV*T1Y(MI)
      T2Y(MI) = ZNV*T1X(MI) - XNV*T1Z(MI)
      T2Z(MI) = XNV*T1Y(MI) - YNV*T1X(MI)
      SALP(MI) = 1.0D0
      IF ( NX.EQ.0 ) GOTO 9
      M = M + NX*NY - 1
      XN2 = X(MI) - S1X - S2X
      YN2 = Y(MI) - S1Y - S2Y
      ZN2 = Z(MI) - S1Z - S2Z
      XXS = T1X(MI)
      YYS = T1Y(MI)
      ZZS = T1Z(MI)
      XT = T2X(MI)
      YT = T2Y(MI)
      ZT = T2Z(MI)
      MI = MI + 1
      DO 8 IY = 1 , NY
         XN2 = XN2 + S2X
         YN2 = YN2 + S2Y
         ZN2 = ZN2 + S2Z
         DO 7 IX = 1 , NX
            XST = IX
            MI = MI - 1
            X(MI) = XN2 + XST*S1X
            Y(MI) = YN2 + XST*S1Y
            Z(MI) = ZN2 + XST*S1Z
            BI(MI) = XA
            SALP(MI) = 1.0D0
            T1X(MI) = XXS
            T1Y(MI) = YYS
            T1Z(MI) = ZZS
            T2X(MI) = XT
            T2Y(MI) = YT
            T2Z(MI) = ZT
    7    CONTINUE
    8 CONTINUE
    9 IPSYM = 0
      NP = N
      MP = M
      RETURN
 
   14 FORMAT (
     &  ' ERROR -- CORNERS OF QUADRILATERAL PATCH DO NOT LIE IN A PLANE'
     &  )
 
      END
