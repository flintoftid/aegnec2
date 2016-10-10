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

      SUBROUTINE PCINT( XI , YI , ZI , CABI , SABI , SALPI , E )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute the interaction matrix elements for the electric field,    **
C* tangent to a segment connected to a surface, due to the current on **
C* the four patches around the connection point.                      **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - XI                                                        **
C* PASSED - YI                                                        **
C* PASSED - ZI                                                        **
C* INPUT  - CABI                                                      **
C* INPUT  - SABI                                                      **
C* INPUT  - SALPI                                                     **
C* OUTPUT - E                                                         **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXK     EXS     S       XJ      YJ      ZJ             **
C* uses value  EXK     EXS     EYK     EYS     EZK     EZS     S      **
C*             T1XJ    T1YJ    T1ZJ    T2XJ    T2YJ    T2ZJ    XJ     **
C*             YJ      ZJ                                             **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       UNERE                                                  **
C* called by   CMSW                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL UNERE

C     Dummy arguments.
      REAL*8 XI , YI , ZI , CABI , SABI , SALPI
      COMPLEX*16 E(9)

C     Local variables.
      INTEGER I1 , I2 , NINT
      REAL*8 DD , DA , DS , F1 , F2 , FCON , G1 , G2 , G3 , G4 , GCON , 
     &       S1 , S2 , S2X , TPI , XXS , XSS , XXJ , XYJ , XZJ , YSS , 
     &       ZSS
      COMPLEX*16 E1 , E2 , E3 , E4 , E5 , E6 , E7 , E8 , E9

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'dataeq.inc'
      
C     Data initialisation.
      DATA TPI /6.283185307179586476925286D0/
      DATA NINT /10/
 

      DD = DSQRT(S)*0.5D0
      DS = 4.0D0*DD/DBLE(NINT)
      DA = DS*DS
      GCON = 1.0D0/S
      FCON = 1.0D0/(2.0D0*TPI*DD)
      XXJ = XJ
      XYJ = YJ
      XZJ = ZJ
      XXS = S
      S = DA
      S1 = DD + DS*0.5D0
      XSS = XJ + S1*(T1XJ+T2XJ)
      YSS = YJ + S1*(T1YJ+T2YJ)
      ZSS = ZJ + S1*(T1ZJ+T2ZJ)
      S1 = S1 + DD
      S2X = S1
      E1 = (0.0D0,0.0D0)
      E2 = (0.0D0,0.0D0)
      E3 = (0.0D0,0.0D0)
      E4 = (0.0D0,0.0D0)
      E5 = (0.0D0,0.0D0)
      E6 = (0.0D0,0.0D0)
      E7 = (0.0D0,0.0D0)
      E8 = (0.0D0,0.0D0)
      E9 = (0.0D0,0.0D0)
      DO 2 I1 = 1 , NINT
         S1 = S1 - DS
         S2 = S2X
         XSS = XSS - DS*T1XJ
         YSS = YSS - DS*T1YJ
         ZSS = ZSS - DS*T1ZJ
         XJ = XSS
         YJ = YSS
         ZJ = ZSS
         DO 1 I2 = 1 , NINT
            S2 = S2 - DS
            XJ = XJ - DS*T2XJ
            YJ = YJ - DS*T2YJ
            ZJ = ZJ - DS*T2ZJ
            CALL UNERE(XI,YI,ZI)
            EXK = EXK*CABI + EYK*SABI + EZK*SALPI
            EXS = EXS*CABI + EYS*SABI + EZS*SALPI
            G1 = (DD+S1)*(DD+S2)*GCON
            G2 = (DD-S1)*(DD+S2)*GCON
            G3 = (DD-S1)*(DD-S2)*GCON
            G4 = (DD+S1)*(DD-S2)*GCON
            F2 = (S1*S1+S2*S2)*TPI
            F1 = S1/F2 - (G1-G2-G3+G4)*FCON
            F2 = S2/F2 - (G1+G2-G3-G4)*FCON
            E1 = E1 + EXK*G1
            E2 = E2 + EXK*G2
            E3 = E3 + EXK*G3
            E4 = E4 + EXK*G4
            E5 = E5 + EXS*G1
            E6 = E6 + EXS*G2
            E7 = E7 + EXS*G3
            E8 = E8 + EXS*G4
            E9 = E9 + EXK*F1 + EXS*F2
    1    CONTINUE
    2 CONTINUE
      E(1) = E1
      E(2) = E2
      E(3) = E3
      E(4) = E4
      E(5) = E5
      E(6) = E6
      E(7) = E7
      E(8) = E8
      E(9) = E9
      XJ = XXJ
      YJ = XYJ
      ZJ = XZJ
      S = XXS
 
      RETURN
 
      END
