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

      SUBROUTINE HANKEL ( Z , H0 , H0P )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Hankel evaluates Hankel function of the first kind, order zero,    **
C* and its derivative for complex argument Z.                         **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - Z                                                         **
C* OUTPUT - H0                                                        **
C* OUTPUT - H0P                                                       **
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
      COMPLEX*16 Z , H0 , H0P

C     Local variables.
      INTEGER I , IB , INIT , IZ , K , MIZ
      INTEGER M(101)
      REAL*8 C1 , C2 , C3 , GAMMA , P0F , P10 , P11 , P20 , P21 , PI , 
     &       PSI , Q10 , Q11 , Q20 , Q21 , TESTT , ZMS
      REAL*8 A1(25), A2(25), A3(25), A4(25), FJX(2)
      COMPLEX*16 CLOGZ , J0 , J0P , P0Z , P1Z , Q0Z , Q1Z , Y0 , Y0P , 
     &           ZI , ZI2 , ZK , FJ

C     Save local variables.
      SAVE M , A1 , A2 , A3 , A4 , FJ , FJX , INIT
C     SAVE PSI , TESTT , ZMS

C     Local equivalencs.
      EQUIVALENCE (FJ,FJX)

C     Data initialisation.
      DATA PI /3.141592654D0/
      DATA GAMMA /0.5772156649D0/
      DATA C1 /-0.0245785095D0/
      DATA C2 /0.3674669052D0/
      DATA C3 /0.7978845608D0/
      DATA P10 /0.0703125D0/
      DATA P20 /0.1121520996D0/
      DATA Q10 /0.125D0/
      DATA Q20 /0.0732421875D0/
      DATA P11 /0.1171875D0/
      DATA P21 /0.1441955566D0/
      DATA Q11 /0.375D0/
      DATA Q21 /0.1025390625D0/
      DATA P0F /0.7853981635D0/
      DATA INIT /0/
      DATA FJX /0.0D0,1.0D0/


      IF (INIT.EQ.0) GO TO 5
1     ZMS=DBLE(Z*DCONJG(Z))
      IF (ZMS.NE.0.0D0) GO TO 2
      PRINT 9
      STOP
2     IB=0
      IF (ZMS.GT.16.81D0) GO TO 4
      IF (ZMS.GT.16.0D0) IB=1

C     Series expansion.

      IZ=INT(1.0D0+ZMS)
      MIZ=M(IZ)
      J0=(1.0D0,0.0D0)
      J0P=J0
      Y0=(0.0D0,0.0D0)
      Y0P=Y0
      ZK=J0
      ZI=Z*Z
      DO 3 K=1,MIZ
      ZK=ZK*A1(K)*ZI
      J0=J0+ZK
      J0P=J0P+A2(K)*ZK
      Y0=Y0+A3(K)*ZK
3     Y0P=Y0P+A4(K)*ZK
      J0P=-0.5D0*Z*J0P
      CLOGZ=CDLOG(0.5D0*Z)
      Y0=(2.0D0*J0*CLOGZ-Y0)/PI+C2
      Y0P=(2.0D0/Z+2.0D0*J0P*CLOGZ+0.5D0*Y0P*Z)/PI+C1*Z
      H0=J0+FJ*Y0
      H0P=J0P+FJ*Y0P
      IF (IB.EQ.0) RETURN
      Y0=H0
      Y0P=H0P

C     Asymptotic expansion.

4     ZI=1.0D0/Z
      ZI2=ZI*ZI
      P0Z=1.0D0+(P20*ZI2-P10)*ZI2
      P1Z=1.0D0+(P11-P21*ZI2)*ZI2
      Q0Z=(Q20*ZI2-Q10)*ZI
      Q1Z=(Q11-Q21*ZI2)*ZI
      ZK=CDEXP(FJ*(Z-P0F))*CDSQRT(ZI)*C3
      H0=ZK*(P0Z+FJ*Q0Z)
      H0P=FJ*ZK*(P1Z+FJ*Q1Z)
      IF (IB.EQ.0) RETURN
      ZMS=DCOS((DSQRT(ZMS)-4.0D0)*31.41592654D0)
      H0=0.5D0*(Y0*(1.0D0+ZMS)+H0*(1.0D0-ZMS))
      H0P=0.5D0*(Y0P*(1.0D0+ZMS)+H0P*(1.0D0-ZMS))
      RETURN

C     Initialization of constants.

5     PSI=-GAMMA
      DO 6 K=1,25
      A1(K)=-0.25D0/(K*K)
      A2(K)=1.0D0/(K+1.0D0)
      PSI=PSI+1.0D0/K
      A3(K)=PSI+PSI
6     A4(K)=(PSI+PSI+1.0D0/(K+1.0D0))/(K+1.0D0)
      DO 8 I=1,101
      TESTT=1.0D0
      DO 7 K=1,24
      INIT=K
      TESTT=-TESTT*I*A1(K)
      IF (TESTT*A3(K).LT.1.0D-6) GO TO 8
7     CONTINUE
8     M(I)=INIT
      GO TO 1
C
9     FORMAT (' ERROR - HANKEL NOT VALID FOR Z=0.')
      END
