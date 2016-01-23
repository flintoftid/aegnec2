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

      SUBROUTINE SUBPH( NX , NY , LD , SALP , M , MP , X , Y , Z , BI , 
     &                  T1X , T1Y , T1Z , T2X , T2Y , T2Z , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* SUBPH generates and modifies patch geometry data.                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NX                                                        **
C* INPUT  - NY                                                        **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
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
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CONECT                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      INTEGER NX , NY , LD , M , MP , DEBUG
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , BI(LD) , T1X(LD) , 
     &       T1Y(LD) , T1Z(LD) , T2X(LD) , T2Y(LD) , T2Z(LD)

C     Local variables.
      INTEGER IX , IY , MI , MIA , NXP , NYP
      REAL*8 S1X , S1Y , S1Z , S2X , S2Y , S2Z , SALN , XA , XXS , 
     &       XST , XT , YYS , YT , ZZS
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SUBPH'
      ENDIF

C     Divide patch for wire connection.
 
      IF ( NY.GT.0 ) GOTO 10
      IF ( NX.EQ.M ) GOTO 10
      NXP = NX + 1
      IX = LD - M
      DO 9 IY = NXP , M
         IX = IX + 1
         NYP = IX - 3
         X(NYP) = X(IX)
         Y(NYP) = Y(IX)
         Z(NYP) = Z(IX)
         BI(NYP) = BI(IX)
         SALP(NYP) = SALP(IX)
         T1X(NYP) = T1X(IX)
         T1Y(NYP) = T1Y(IX)
         T1Z(NYP) = T1Z(IX)
         T2X(NYP) = T2X(IX)
         T2Y(NYP) = T2Y(IX)
         T2Z(NYP) = T2Z(IX)
    9 CONTINUE
   10 MI = LD + 1 - NX
      XXS = X(MI)
      YYS = Y(MI)
      ZZS = Z(MI)
      XA = BI(MI)*0.25D0
      XST = DSQRT(XA)*0.5D0
      S1X = T1X(MI)
      S1Y = T1Y(MI)
      S1Z = T1Z(MI)
      S2X = T2X(MI)
      S2Y = T2Y(MI)
      S2Z = T2Z(MI)
      SALN = SALP(MI)
      XT = XST
      YT = XST
      IF ( NY.GT.0 ) GOTO 11
      MIA = MI
      GOTO 12
   11 M = M + 1
      MP = MP + 1
      MIA = LD + 1 - M
   12 DO 13 IX = 1 , 4
         X(MIA) = XXS + XT*S1X + YT*S2X
         Y(MIA) = YYS + XT*S1Y + YT*S2Y
         Z(MIA) = ZZS + XT*S1Z + YT*S2Z
         BI(MIA) = XA
         T1X(MIA) = S1X
         T1Y(MIA) = S1Y
         T1Z(MIA) = S1Z
         T2X(MIA) = S2X
         T2Y(MIA) = S2Y
         T2Z(MIA) = S2Z
         SALP(MIA) = SALN
         IF ( IX.EQ.2 ) YT = -YT
         IF ( IX.EQ.1 .OR. IX.EQ.3 ) XT = -XT
         MIA = MIA - 1
   13 CONTINUE
      M = M + 3
      IF ( NX.LE.MP ) MP = MP + 3
      IF ( NY.GT.0 ) Z(MI) = 10000.0D0
 
      RETURN
 
C  14 FORMAT (
C    &  ' ERROR -- CORNERS OF QUADRILATERAL PATCH DO NOT LIE IN A PLANE'
C    &  )
 
      END
