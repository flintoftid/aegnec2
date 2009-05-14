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

      SUBROUTINE REFLC( IX , IY , IZ , ITX , NOP , LD , SALP ,
     &                  IPSYM , N1 , N2 , N , NP , M1 , M2 , M , MP , 
     &                  ITAG , X , Y , Z , BI , X2 , Y2 , Z2 , T1X , 
     &                  T1Y , T1Z , T2X , T2Y , T2Z , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* REFLC reflects partial structure along X,Y, or Z axes or rotates   **
C* structure to complete a symmetric structure.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - IX                                                        **
C* INPUT  - IY                                                        **
C* INPUT  - IZ                                                        **
C* INPUT  - ITX                                                       **
C* INPUT  - NOP                                                       **
C* INPUT  - LD                                                        **
C* OUTPUT - SALP                                                      **
C* OUTPUT - IPSYM                                                     **
C* INPUT  - N1                                                        **
C* INPUT  - N2                                                        **
C* OUTPUT - N                                                         **
C* OUTPUT - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M2                                                        **
C* OUTPUT - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - ITAG                                                      **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* OUTPUT - BI                                                        **
C* OUTPUT - X2                                                        **
C* OUTPUT - Y2                                                        **
C* OUTPUT - Z2                                                        **
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

C     Dummy arguments.
      INTEGER IX , IY , IZ , ITX , NOP , LD , IPSYM , N1 , N2 , N , NP , 
     &        M1 , M2 , M , MP , IFAIL , DEBUG
      INTEGER ITAG(2*LD)
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , BI(LD) , X2(LD) , 
     &       Y2(LD) , Z2(LD) , T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD)

C     Local variables.
      INTEGER I , ITAGI , ITI , J , K , NX , NXX
      REAL*8 CS , E1 , E2 , FNOP , SAM , SS , XK , YK


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING REFLC'
      ENDIF
 
      NP = N
      MP = M
      IPSYM = 0
      ITI = ITX
      IF ( IX.LT.0 ) GOTO 19
      IF ( NOP.EQ.0 ) RETURN
      IPSYM = 1
      IF ( IZ.EQ.0 ) GOTO 6

C     Reflect along z axis.

      IPSYM = 2
      IF ( N.LT.N2 ) GOTO 3
      DO 2 I = N2 , N
         NX = I + N - N1
         E1 = Z(I)
         E2 = Z2(I)
         IF ( DABS(E1)+DABS(E2).GT.1.D-5 .AND. E1*E2.GE.-1.D-6 ) GOTO 1
         WRITE (3,24) I
         IFAIL=30
         RETURN
    1    X(NX) = X(I)
         Y(NX) = Y(I)
         Z(NX) = -E1
         X2(NX) = X2(I)
         Y2(NX) = Y2(I)
         Z2(NX) = -E2
         ITAGI = ITAG(I)
         IF ( ITAGI.EQ.0 ) ITAG(NX) = 0
         IF ( ITAGI.NE.0 ) ITAG(NX) = ITAGI + ITI
         BI(NX) = BI(I)
    2 CONTINUE
      N = N*2 - N1
      ITI = ITI*2
    3 IF ( M.LT.M2 ) GOTO 6
      NXX = LD + 1 - M1
      DO 5 I = M2 , M
         NXX = NXX - 1
         NX = NXX - M + M1
         IF ( DABS(Z(NXX)).GT.1.D-10 ) GOTO 4
         WRITE (3,25) I
         IFAIL=31
         RETURN
    4    X(NX) = X(NXX)
         Y(NX) = Y(NXX)
         Z(NX) = -Z(NXX)
         T1X(NX) = T1X(NXX)
         T1Y(NX) = T1Y(NXX)
         T1Z(NX) = -T1Z(NXX)
         T2X(NX) = T2X(NXX)
         T2Y(NX) = T2Y(NXX)
         T2Z(NX) = -T2Z(NXX)
         SALP(NX) = -SALP(NXX)
         BI(NX) = BI(NXX)
    5 CONTINUE
      M = M*2 - M1
    6 IF ( IY.EQ.0 ) GOTO 12

C     Reflect along y axis.

      IF ( N.LT.N2 ) GOTO 9
      DO 8 I = N2 , N
         NX = I + N - N1
         E1 = Y(I)
         E2 = Y2(I)
         IF ( DABS(E1)+DABS(E2).GT.1.D-5 .AND. E1*E2.GE.-1.D-6 ) GOTO 7
         WRITE (3,24) I
         IFAIL=30
         RETURN
    7    X(NX) = X(I)
         Y(NX) = -E1
         Z(NX) = Z(I)
         X2(NX) = X2(I)
         Y2(NX) = -E2
         Z2(NX) = Z2(I)
         ITAGI = ITAG(I)
         IF ( ITAGI.EQ.0 ) ITAG(NX) = 0
         IF ( ITAGI.NE.0 ) ITAG(NX) = ITAGI + ITI
         BI(NX) = BI(I)
    8 CONTINUE
      N = N*2 - N1
      ITI = ITI*2
    9 IF ( M.LT.M2 ) GOTO 12
      NXX = LD + 1 - M1
      DO 11 I = M2 , M
         NXX = NXX - 1
         NX = NXX - M + M1
         IF ( DABS(Y(NXX)).GT.1.D-10 ) GOTO 10
         WRITE (3,25) I
         IFAIL=31
         RETURN
   10    X(NX) = X(NXX)
         Y(NX) = -Y(NXX)
         Z(NX) = Z(NXX)
         T1X(NX) = T1X(NXX)
         T1Y(NX) = -T1Y(NXX)
         T1Z(NX) = T1Z(NXX)
         T2X(NX) = T2X(NXX)
         T2Y(NX) = -T2Y(NXX)
         T2Z(NX) = T2Z(NXX)
         SALP(NX) = -SALP(NXX)
         BI(NX) = BI(NXX)
   11 CONTINUE
      M = M*2 - M1
   12 IF ( IX.EQ.0 ) GOTO 18

C     Reflect along x axis.

      IF ( N.LT.N2 ) GOTO 15
      DO 14 I = N2 , N
         NX = I + N - N1
         E1 = X(I)
         E2 = X2(I)
         IF ( DABS(E1)+DABS(E2).GT.1.D-5 .AND. E1*E2.GE.-1.D-6 ) GOTO 13
         WRITE (3,24) I
         IFAIL=30
         RETURN
   13    X(NX) = -E1
         Y(NX) = Y(I)
         Z(NX) = Z(I)
         X2(NX) = -E2
         Y2(NX) = Y2(I)
         Z2(NX) = Z2(I)
         ITAGI = ITAG(I)
         IF ( ITAGI.EQ.0 ) ITAG(NX) = 0
         IF ( ITAGI.NE.0 ) ITAG(NX) = ITAGI + ITI
         BI(NX) = BI(I)
   14 CONTINUE
      N = N*2 - N1
   15 IF ( M.LT.M2 ) GOTO 18
      NXX = LD + 1 - M1
      DO 17 I = M2 , M
         NXX = NXX - 1
         NX = NXX - M + M1
         IF ( DABS(X(NXX)).GT.1.D-10 ) GOTO 16
         WRITE (3,25) I
         IFAIL=31
         RETURN
   16    X(NX) = -X(NXX)
         Y(NX) = Y(NXX)
         Z(NX) = Z(NXX)
         T1X(NX) = -T1X(NXX)
         T1Y(NX) = T1Y(NXX)
         T1Z(NX) = T1Z(NXX)
         T2X(NX) = -T2X(NXX)
         T2Y(NX) = T2Y(NXX)
         T2Z(NX) = T2Z(NXX)
         SALP(NX) = -SALP(NXX)
         BI(NX) = BI(NXX)
   17 CONTINUE
      M = M*2 - M1
   18 RETURN

C     Reproduce structure with rotation to form cylindrical structure.

   19 FNOP = NOP
      IPSYM = -1
      SAM = 6.283185307179586476925286D0/FNOP
      CS = DCOS(SAM)
      SS = DSIN(SAM)
      IF ( N.LT.N2 ) GOTO 21
      N = N1 + (N-N1)*NOP
      NX = NP + 1
      DO 20 I = NX , N
         K = I - NP + N1
         XK = X(K)
         YK = Y(K)
         X(I) = XK*CS - YK*SS
         Y(I) = XK*SS + YK*CS
         Z(I) = Z(K)
         XK = X2(K)
         YK = Y2(K)
         X2(I) = XK*CS - YK*SS
         Y2(I) = XK*SS + YK*CS
         Z2(I) = Z2(K)
         ITAGI = ITAG(K)
         IF ( ITAGI.EQ.0 ) ITAG(I) = 0
         IF ( ITAGI.NE.0 ) ITAG(I) = ITAGI + ITI
         BI(I) = BI(K)
   20 CONTINUE
   21 IF ( M.LT.M2 ) GOTO 23
      M = M1 + (M-M1)*NOP
      NX = MP + 1
      K = LD + 1 - M1
      DO 22 I = NX , M
         K = K - 1
         J = K - MP + M1
         XK = X(K)
         YK = Y(K)
         X(J) = XK*CS - YK*SS
         Y(J) = XK*SS + YK*CS
         Z(J) = Z(K)
         XK = T1X(K)
         YK = T1Y(K)
         T1X(J) = XK*CS - YK*SS
         T1Y(J) = XK*SS + YK*CS
         T1Z(J) = T1Z(K)
         XK = T2X(K)
         YK = T2Y(K)
         T2X(J) = XK*CS - YK*SS
         T2Y(J) = XK*SS + YK*CS
         T2Z(J) = T2Z(K)
         SALP(J) = SALP(K)
         BI(J) = BI(K)
   22 CONTINUE
 
   23 RETURN
 
   24 FORMAT (' GEOMETRY DATA ERROR--SEGMENT',I5,
     &        ' LIES IN PLANE OF SYMMETRY')
   25 FORMAT (' GEOMETRY DATA ERROR--PATCH',I4,
     &        ' LIES IN PLANE OF SYMMETRY')
 
      END

