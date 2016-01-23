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
 
      SUBROUTINE HELIX( SSS , HL , A1 , B1 , A2 , B2 , RAD , NS , ITG , 
     &                  LD , IPSYM , N , NP , M , MP , ITAG , X , 
     &                  Y , Z , BI , X2 , Y2 , Z2 , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Subroutine HELIX generates segment geometry data for a helix of NS **
C* segments.                                                          **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - SSS                                                       **
C* INPUT  - HL                                                        **
C* INPUT  - A1                                                        **
C* OUTPUT - B1                                                        **
C* INPUT  - A2                                                        **
C* OUTPUT - B2                                                        **
C* INPUT  - RAD                                                       **
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

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NS , ITG , LD , IPSYM , N , NP , M , MP , DEBUG 
      INTEGER ITAG(2*LD)
      REAL*8 SSS , HL , A1 , A2 , B1 , B2 , RAD 
      REAL*8 X(LD) , Y(LD) , Z(LD) , BI(LD) , X2(LD) , Y2(LD) , Z2(LD)

C     Local variables.
      INTEGER I , IST
      REAL*8 COPY , HDIA , HMAJ , HMIN , PI , PITCH , SANGLE , TURN , 
     &       ZINC
 
C     Data initialisation.
      DATA PI /3.141592653589793238462643D0/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING HELIX'
      ENDIF
 
      IST = N + 1
      N = N + NS
      NP = N
      MP = M
      IPSYM = 0
      IF ( NS.LT.1 ) RETURN
      ZINC = DABS(HL/NS)
      Z(IST) = 0.0D0
      DO 25 I = IST , N
         BI(I) = RAD
         ITAG(I) = ITG
         IF ( I.NE.IST ) Z(I) = Z(I-1) + ZINC
         Z2(I) = Z(I) + ZINC
         IF ( A2.NE.A1 ) GOTO 10
         IF ( B1.EQ.0.0D0 ) B1 = A1
         X(I) = A1*DCOS(2.0D0*PI*Z(I)/SSS)
         Y(I) = B1*DSIN(2.0D0*PI*Z(I)/SSS)
         X2(I) = A1*DCOS(2.0D0*PI*Z2(I)/SSS)
         Y2(I) = B1*DSIN(2.0D0*PI*Z2(I)/SSS)
         GOTO 20
   10    IF ( B2.EQ.0.0D0 ) B2 = A2
         X(I) = (A1+(A2-A1)*Z(I)/DABS(HL))*DCOS(2.0D0*PI*Z(I)/SSS)
         Y(I) = (B1+(B2-B1)*Z(I)/DABS(HL))*DSIN(2.0D0*PI*Z(I)/SSS)
         X2(I) = (A1+(A2-A1)*Z2(I)/DABS(HL))*DCOS(2.0D0*PI*Z2(I)/SSS)
         Y2(I) = (B1+(B2-B1)*Z2(I)/DABS(HL))*DSIN(2.0D0*PI*Z2(I)/SSS)
   20    IF ( HL.GT.0 ) GOTO 25
         COPY = X(I)
         X(I) = Y(I)
         Y(I) = COPY
         COPY = X2(I)
         X2(I) = Y2(I)
         Y2(I) = COPY
   25 CONTINUE
      IF ( A2.EQ.A1 ) GOTO 21
      SANGLE = DATAN(A2/(DABS(HL)+(DABS(HL)*A1)/(A2-A1)))
      WRITE (CHRSLT,104) SANGLE
      RETURN
   21 IF ( A1.NE.B1 ) GOTO 30
      HDIA = 2.0D0*A1
      TURN = HDIA*PI
      PITCH = DATAN(SSS/(PI*HDIA))
      TURN = TURN/DCOS(PITCH)
      PITCH = 180.0D0*PITCH/PI
      GOTO 40
   30 IF ( A1.LT.B1 ) GOTO 34
      HMAJ = 2.0D0*A1
      HMIN = 2.0D0*B1
      GOTO 35
   34 HMAJ = 2.0D0*B1
      HMIN = 2.0D0*A1
   35 HDIA = DSQRT((HMAJ**2+HMIN**2)/2*HMAJ)
      TURN = 2.0D0*PI*HDIA
      PITCH = (180.0D0/PI)*DATAN(SSS/(PI*HDIA))
   40 WRITE (CHRSLT,105) PITCH , TURN
 
      RETURN
  104 FORMAT (5X,'THE CONE ANGLE OF THE SPIRAL IS',F10.4)
  105 FORMAT (5X,'THE PITCH ANGLE IS',F10.4/5X,
     &        'THE LENGTH OF WIRE/TURN IS',F10.4)
 
      END
