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

      COMPLEX*16 FUNCTION FBAR( P )

C*--------------------------------------------------------------------**
C*                                                                    **
C* FBAR is Sommerfeld attenuation function for Norton's asymptotic    **
C* field approximation at numerical distance P.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - P                                                         **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   GWAVE                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Dummy arguments.
      COMPLEX*16 P

C     Local variables.
      INTEGER I , MINUS
      REAL*8 ACCS , SMS , SP , TMS , TOSP
      REAL*8 FJX(2)
      COMPLEX*16 ZZ , ZZS , SUM , POW , TERM , FJ

C     Equivalent storage.
      EQUIVALENCE (FJ,FJX)
 
C     Data initialisation.
      DATA TOSP /1.128379167095512573896158D0/
      DATA ACCS /1.0D-12/
      DATA SP /1.772453850905516027298167D0/
      DATA FJX/0.0D0 , 1.0D0/
 

      ZZ = FJ*CDSQRT(P)
      IF ( CDABS(ZZ).GT.3.0D0 ) GOTO 3

C     Series expansion.

      ZZS = ZZ*ZZ
      SUM = ZZ
      POW = ZZ
      DO 1 I = 1 , 100
         POW = -POW*ZZS/DBLE(I)
         TERM = POW/(2.0D0*I+1.0D0)
         SUM = SUM + TERM
         TMS = DBLE(TERM*DCONJG(TERM))
         SMS = DBLE(SUM*DCONJG(SUM))
         IF ( TMS/SMS.LT.ACCS ) GOTO 2
    1 CONTINUE
    2 FBAR = 1.0D0 - (1.0D0-SUM*TOSP)*ZZ*CDEXP(ZZS)*SP
      RETURN

C     Asymptotic expansion.

    3 IF ( DBLE(ZZ).GE.0. ) GOTO 4
      MINUS = 1
      ZZ = -ZZ
      GOTO 5
    4 MINUS = 0
    5 ZZS = 0.5D0/(ZZ*ZZ)
      SUM = DCMPLX(0.0D0,0.0D0)
      TERM = DCMPLX(1.0D0,0.0D0)
      DO 6 I = 1 , 6
         TERM = -TERM*(2.0D0*I-1.0D0)*ZZS
         SUM = SUM + TERM
    6 CONTINUE
      IF ( MINUS.EQ.1 ) SUM = SUM - 2.0D0*SP*ZZ*CDEXP(ZZ*ZZ)
      FBAR = -SUM
 
      RETURN
 
      END
