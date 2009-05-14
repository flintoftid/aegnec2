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

      SUBROUTINE BESSEL ( Z , J0 , J0P )

C*--------------------------------------------------------------------**
C*                                                                    **
C* BESSEL evaluates the zero-order bessel function and its derivative **
C* for complex argument Z.                                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - Z                                                         **
C* OUTPUT - J0                                                        **
C* OUTPUT - J0P                                                       **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   SAOA                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Dummy arguments.
      COMPLEX*16 Z , J0 , J0P

C     Local variables.
      INTEGER I , IB , INIT , IZ , K , MIZ
      INTEGER M(101)
      REAL*8 C3 , P10 , P11 , P20 , P21 , POF , Q10 , Q11 , Q20 , 
     &       Q21 , TESTT , ZMS
      REAL*8 A1(25) , A2(25) , FJX(2)
      COMPLEX*16 P0Z , P1Z , Q0Z , Q1Z , ZI , ZI2 , ZK , FJ , CZ , 
     &           SZ , J0X , J0PX

C     Save local variables.
      SAVE M , A1 , A2 , FJX , FJ , INIT

C     Local equivalences.
      EQUIVALENCE (FJ,FJX)

C     Data initialisation.
      DATA C3  /0.7978845608D0/
      DATA P10 /0.0703125D0/
      DATA P20 /0.1121520996D0/
      DATA Q10 /0.125D0/
      DATA Q20 /0.0732421875D0/
      DATA P11 /0.1171875D0/
      DATA P21 /0.1441955566D0/
      DATA Q11 /0.375D0/
      DATA Q21 /0.1025390625D0/
      DATA POF  /0.7853981635D0/
      DATA INIT /0/
      DATA FJX  /0.0D0,1.0D0/


      IF (INIT.EQ.0) GO TO 5
1     ZMS=DBLE(Z*DCONJG(Z))
      IF (ZMS.GT.1.0D-12) GO TO 2
      J0=(1.0D0,0.0D0)
      J0P=-0.5D0*Z
      RETURN
2     IB=0
      IF (ZMS.GT.37.21D0) GO TO 4
      IF (ZMS.GT.36.0D0) IB=1

C     Series expansion.
      IZ=INT(1.0D0+ZMS)
      MIZ=M(IZ)
      J0=(1.0D0,0.0D0)
      J0P=J0
      ZK=J0
      ZI=Z*Z
      DO 3 K=1,MIZ
      ZK=ZK*A1(K)*ZI
      J0=J0+ZK
3     J0P=J0P+A2(K)*ZK
      J0P=-0.5D0*Z*J0P
      IF (IB.EQ.0) RETURN
      J0X=J0
      J0PX=J0P

C     Asymptotic expansion.

4     ZI=1.0D0/Z
      ZI2=ZI*ZI
      P0Z=1.0D0+(P20*ZI2-P10)*ZI2
      P1Z=1.0D0+(P11-P21*ZI2)*ZI2
      Q0Z=(Q20*ZI2-Q10)*ZI
      Q1Z=(Q11-Q21*ZI2)*ZI
      ZK=CDEXP(FJ*(Z-POF))
      ZI2=1.0D0/ZK
      CZ=0.5D0*(ZK+ZI2)
      SZ=FJ*0.5D0*(ZI2-ZK)
      ZK=C3*CDSQRT(ZI)
      J0=ZK*(P0Z*CZ-Q0Z*SZ)
      J0P=-ZK*(P1Z*SZ+Q1Z*CZ)
      IF (IB.EQ.0) RETURN
      ZMS=DCOS((SQRT(ZMS)-6.0D0)*31.41592654D0)
      J0=0.5D0*(J0X*(1.0D0+ZMS)+J0*(1.0D0-ZMS))
      J0P=0.5D0*(J0PX*(1.0D0+ZMS)+J0P*(1.0D0-ZMS))
      RETURN

C     Initialization of constants.

5     DO 6 K=1,25
      A1(K)=-0.25D0/(K*K)
6     A2(K)=1.0D0/(K+1.0D0)
      DO 8 I=1,101
      TESTT=1.0D0
      DO 7 K=1,24
      INIT=K
      TESTT=-TESTT*I*A1(K)
      IF (TESTT.LT.1.0D-6) GO TO 8
7     CONTINUE
8     M(I)=INIT
      GO TO 1

      END
