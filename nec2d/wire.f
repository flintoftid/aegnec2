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

      SUBROUTINE WIRE( XW1 , YW1 , ZW1 , XW2 , YW2 , ZW2 , RAD ,
     &                 RDEL , RRAD , NS , ITG , LD , IPSYM , 
     &                 N , NP , M , MP , ITAG , X , Y , Z , BI , X2 , 
     &                 Y2 , Z2 , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine WIRE generates segment geometry data for a straight     **
C* wire of NS segments.                                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XW1                                                       **
C* INPUT  - YW1                                                       **
C* INPUT  - ZW1                                                       **
C* INPUT  - XW2                                                       **
C* INPUT  - YW2                                                       **
C* INPUT  - ZW2                                                       **
C* INPUT  - RAD                                                       **
C* INPUT  - RDEL                                                      **
C* INPUT  - RRAD                                                      **
C* INPUT  - NS                                                        **
C* INPUT  - ITG                                                       **
C* INPUT  - LD                                                        **
C* OUTPUT - IPSYM                                                     **
C* OUTPUT - N                                                         **
C* OUTPUT - NP                                                        **
C* INPUT  - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - ITAG                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - BI                                                        **
C* OUTPUT - X2                                                        **
C* OUTPUT - Y2                                                        **
C* OUTPUT - Z2                                                        **
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

C     Dummy arguments.
      INTEGER NS , ITG , LD , IPSYM , N , NP , M , MP , DEBUG
      INTEGER ITAG(2*LD)
      REAL*8 XW1 , YW1 , ZW1 , XW2 , YW2 , ZW2 , RAD , RDEL , RRAD
      REAL*8 X(LD) , Y(LD) , Z(LD) , BI(LD) , 
     &       X2(LD) , Y2(LD) , Z2(LD)

C     Local variables.
      INTEGER I , IST
      REAL*8 DELZ , FNS , RADZ , RD , XD , XS1 , XXS2 , YD , YS1 , 
     &       YS2 , ZD , ZS1 , ZS2
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING WIRE'
      ENDIF
 
      IST = N + 1
      N = N + NS
      NP = N
      MP = M
      IPSYM = 0
      IF ( NS.LT.1 ) RETURN
      XD = XW2 - XW1
      YD = YW2 - YW1
      ZD = ZW2 - ZW1
      IF ( DABS(RDEL-1.0D0).LT.1.D-6 ) GOTO 1
      DELZ = DSQRT(XD*XD+YD*YD+ZD*ZD)
      XD = XD/DELZ
      YD = YD/DELZ
      ZD = ZD/DELZ
      DELZ = DELZ*(1.0D0-RDEL)/(1.0D0-RDEL**NS)
      RD = RDEL
      GOTO 2
    1 FNS = NS
      XD = XD/FNS
      YD = YD/FNS
      ZD = ZD/FNS
      DELZ = 1.0D0
      RD = 1.0D0
    2 RADZ = RAD
      XS1 = XW1
      YS1 = YW1
      ZS1 = ZW1
      DO 3 I = IST , N
         ITAG(I) = ITG
         XXS2 = XS1 + XD*DELZ
         YS2 = YS1 + YD*DELZ
         ZS2 = ZS1 + ZD*DELZ
         X(I) = XS1
         Y(I) = YS1
         Z(I) = ZS1
         X2(I) = XXS2
         Y2(I) = YS2
         Z2(I) = ZS2
         BI(I) = RADZ
         DELZ = DELZ*RD
         RADZ = RADZ*RRAD
         XS1 = XXS2
         YS1 = YS2
         ZS1 = ZS2
    3 CONTINUE
      X2(N) = XW2
      Y2(N) = YW2
      Z2(N) = ZW2
 
      RETURN
 
      END
