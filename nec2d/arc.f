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

      SUBROUTINE ARC ( ITG , NS , RADA , ANG1 , ANG2 , RAD , LD , 
     &                 IPSYM , N , NP , M , MP , ITAG , X , Y , Z , BI , 
     &                 X2 , Y2 , Z2 , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* ARC fills arrays with segment corrdinates for a circular arc of    **
C* segments.                                                          **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ITG                                                       **
C* INPUT  - NS                                                        **
C* INPUT  - RADA                                                      **
C* INPUT  - ANG1                                                      **
C* INPUT  - ANG2                                                      **
C* INPUT  - RAD                                                       **
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
      INTEGER ITG , NS , LD , IPSYM , N , NP , M , MP , IFAIL , DEBUG
      INTEGER ITAG(2*LD)
      REAL*8 RADA , ANG1 , ANG2 , RAD
      REAL*8 X(LD) , Y(LD) , Z(LD) , BI(LD) , X2(LD) , Y2(LD) , Z2(LD)

C     Local variables.
      REAL*8 ANG , DANG , TA , XS1 , XXS2 , ZS1 , ZS2, ANGTOL
      INTEGER I , IST 

C     Data initialisation.
      DATA TA/0.0174532925199432957692369D0/
      DATA ANGTOL/360.00000001D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING ARC'
      ENDIF

      IST = N + 1
      N = N + NS
      NP = N
      MP = M
      IPSYM = 0
      IF ( NS.LT.1 ) RETURN
      IF ( DABS(ANG2-ANG1).LT.ANGTOL ) GOTO 1
      WRITE (CHRSLT,3)
      IFAIL=2
      RETURN
    1 ANG = ANG1*TA
      DANG = (ANG2-ANG1)*TA/NS
      XS1 = RADA*DCOS(ANG)
      ZS1 = RADA*DSIN(ANG)
      DO 2 I = IST , N
         ANG = ANG + DANG
         XXS2 = RADA*DCOS(ANG)
         ZS2 = RADA*DSIN(ANG)
         X(I) = XS1
         Y(I) = 0.0D0
         Z(I) = ZS1
         X2(I) = XXS2
         Y2(I) = 0.0D0
         Z2(I) = ZS2
         XS1 = XXS2
         ZS1 = ZS2
         BI(I) = RAD
         ITAG(I) = ITG
    2 CONTINUE

      RETURN

    3 FORMAT (' ERROR -- ARC ANGLE EXCEEDS 360. DEGREES')

      END
