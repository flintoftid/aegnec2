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

      SUBROUTINE CMSS( J1 , J2 , IM1 , IM2 , CCM , NROW , ITRP , 
     &                 LD , SALP , X , Y , Z , BI , T1X , T1Y , T1Z , 
     &                 T2X , T2Y , T2Z , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CMSS computes and stores matrix elements for the H field at patch  **
C* centers due to the current on patches.                             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - J1                                                        **
C* INPUT  - J2                                                        **
C* INPUT  - IM1                                                       **
C* INPUT  - IM2                                                       **
C* OUTPUT - CCM                                                       **
C* INPUT  - NROW                                                      **
C* INPUT  - ITRP                                                      **
C* INPUT  - LD                                                        **
C* INPUT  - SALP                                                      **
C* INPUT  - X                                                         **
C* INPUT  - Y                                                         **
C* INPUT  - Z                                                         **
C* INPUT  - BI                                                        **
C* INPUT  - T1X                                                       **
C* INPUT  - T1Y                                                       **
C* INPUT  - T1Z                                                       **
C* INPUT  - T2X                                                       **
C* INPUT  - T2Y                                                       **
C* INPUT  - T2Z                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    S       T1XJ    T1YJ    T1ZJ    T2XJ    T2YJ    T2ZJ   **
C*             XJ      YJ      ZJ                                     **
C* uses value  EXK     EXS     EYK     EYS     EZK     EZS            **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       HINTG                                                  **
C* called by   CMNGF   CMSET                                          **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL HINTG

C     Dummy arguments.
      INTEGER J1 , J2 , IM1 , IM2 , NROW , ITRP , LD , DEBUG
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       BI(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , 
     &       T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 CCM(NROW,1)

C     Local variables.
      INTEGER I , I1 , I2 , ICOMP , II1 , II2 , IL , J , JJ1 , JJ2 , 
     &        JL , LDP 
      REAL*8 T1XI , T1YI , T1ZI , T2XI , T2YI , T2ZI , XI , YI , ZI
      COMPLEX*16 G11 , G12 , G21 , G22

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMSS'
      ENDIF
 
      LDP = LD + 1
      I1 = (IM1+1)/2
      I2 = (IM2+1)/2
      ICOMP = I1*2 - 3
      II1 = -1
      IF ( ICOMP+2.LT.IM1 ) II1 = -2
C     Loop over observation patches.
      DO 6 I = I1 , I2
         IL = LDP - I
         ICOMP = ICOMP + 2
         II1 = II1 + 2
         II2 = II1 + 1
         T1XI = T1X(IL)*SALP(IL)
         T1YI = T1Y(IL)*SALP(IL)
         T1ZI = T1Z(IL)*SALP(IL)
         T2XI = T2X(IL)*SALP(IL)
         T2YI = T2Y(IL)*SALP(IL)
         T2ZI = T2Z(IL)*SALP(IL)
         XI = X(IL)
         YI = Y(IL)
         ZI = Z(IL)
         JJ1 = -1
C     Loop over source patches.
         DO 5 J = J1 , J2
            JL = LDP - J
            JJ1 = JJ1 + 2
            JJ2 = JJ1 + 1
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
            CALL HINTG(XI,YI,ZI,DEBUG)
            G11 = -(T2XI*EXK+T2YI*EYK+T2ZI*EZK)
            G12 = -(T2XI*EXS+T2YI*EYS+T2ZI*EZS)
            G21 = -(T1XI*EXK+T1YI*EYK+T1ZI*EZK)
            G22 = -(T1XI*EXS+T1YI*EYS+T1ZI*EZS)
            IF ( I.NE.J ) GOTO 1
            G11 = G11 - 0.5D0
            G22 = G22 + 0.5D0
    1       IF ( ITRP.NE.0 ) GOTO 3
C     Normal fill.
            IF ( ICOMP.LT.IM1 ) GOTO 2
            CCM(II1,JJ1) = G11
            CCM(II1,JJ2) = G12
    2       IF ( ICOMP.GE.IM2 ) GOTO 5
            CCM(II2,JJ1) = G21
            CCM(II2,JJ2) = G22
            GOTO 5
C     Transposed fill.
    3       IF ( ICOMP.LT.IM1 ) GOTO 4
            CCM(JJ1,II1) = G11
            CCM(JJ2,II1) = G12
    4       IF ( ICOMP.GE.IM2 ) GOTO 5
            CCM(JJ1,II2) = G21
            CCM(JJ2,II2) = G22
    5    CONTINUE
    6 CONTINUE
 
      RETURN
 
      END
