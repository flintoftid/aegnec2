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
 
      SUBROUTINE CMSET( NROW , CCM , RKHX , IEXKX , LD , SALP ,
     &                  NLOAD , ZARRAY , D , N1 , N , NP , M1 , M , MP , 
     &                  ICON1 , ICON2 , ICONX , X , Y , Z , SI , BI , 
     &                  CAB , SAB , T1X , T1Y , T1Z , T2X , T2Y , T2Z , 
     &                  JMAX , JSNO , JCO , AX , BX , CX , IFAIL , 
     &                  DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CMSET sets up the complex structure matrix in the array CCM.       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - NROW                                                      **
C* OUTPUT - CCM                                                       **
C* INPUT  - RKHX                                                      **
C* INPUT  - IEXKX                                                     **
C* INPUT  - LD                                                        **
C* PASSED - SALP                                                      **
C* INPUT  - NLOAD                                                     **
C* INPUT  - ZARRAY                                                    **
C* OUTPUT - D                                                         **
C* PASSED - N1                                                        **
C* INPUT  - N                                                         **
C* INPUT  - NP                                                        **
C* PASSED - M1                                                        **
C* INPUT  - M                                                         **
C* INPUT  - MP                                                        **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* PASSED - ICONX                                                     **
C* PASSED - X                                                         **
C* PASSED - Y                                                         **
C* PASSED - Z                                                         **
C* PASSED - SI                                                        **
C* PASSED - BI                                                        **
C* PASSED - CAB                                                       **
C* PASSED - SAB                                                       **
C* PASSED - T1X                                                       **
C* PASSED - T1Y                                                       **
C* PASSED - T1Z                                                       **
C* PASSED - T2X                                                       **
C* PASSED - T2Y                                                       **
C* PASSED - T2Z                                                       **
C* INPUT  - JMAX                                                      **
C* INPUT  - JSNO                                                      **
C* OUTPUT - JCO                                                       **
C* INPUT  - AX                                                        **
C* PASSED - BX                                                        **
C* INPUT  - CX                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    IEXK    RKH                                            **
C* uses value  ICASE   NBLOKS  NLAST   NPBLK   SSX                    **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       BLCKOT  CMSS    CMSW    CMWS    CMWW    TRIO           **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL BLCKOT , CMSS , CMSW , CMWS , CMWW , TRIO 

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NROW , IEXKX , LD , NLOAD , N1 , N , NP , M1 , M , MP , 
     &        JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , ICONX(LD) , 
     &        JCO(JMAX)
      REAL*8 RKHX
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) , BI(LD) , CAB(LD) , SAB(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX)  
      COMPLEX*16 CCM(NROW,1) , ZARRAY(LD) , D(2*LD)

C     Local variables.
      INTEGER I , I1 , I2 , IJ , IM1 , IM2 , IN2 , IOUT , IPR , IST ,
     &        ISV , IT , IXBLK1 , J , JM1 , JM2 , JSS , JST , K , 
     &        KA , KK , MP2 , NNEQ , NOP , NNPEQ 
      COMPLEX*16 ZAJ , DETER

C     Common storage.
      INCLUDE 'matpar.inc'
      INCLUDE 'smat.inc'
      INCLUDE 'dataj.inc'
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMSET'
      ENDIF
 
      MP2 = 2*MP
      NNPEQ = NP + MP2
      NNEQ = N + 2*M
      NOP = NNEQ/NNPEQ
      IF ( ICASE.GT.2 ) REWIND CHTMP1
C     In NEC81 following two lines are commented out.
      RKH = RKHX
      IEXK = IEXKX
C
      IOUT = 2*NPBLK*NROW
      IT = NPBLK

C     Cycle over matrix blocks.

      DO 16 IXBLK1 = 1 , NBLOKS
         ISV = (IXBLK1-1)*NPBLK
         IF ( IXBLK1.EQ.NBLOKS ) IT = NLAST
         DO 2 I = 1 , NROW
            DO 1 J = 1 , IT
               CCM(I,J) = (0.0D0,0.0D0)
    1       CONTINUE
    2    CONTINUE
         I1 = ISV + 1
         I2 = ISV + IT
         IN2 = I2
         IF ( IN2.GT.NP ) IN2 = NP
         IM1 = I1 - NP
         IM2 = I2 - NP
         IF ( IM1.LT.1 ) IM1 = 1
         IST = 1
         IF ( I1.LE.NP ) IST = NP - I1 + 2
         IF ( N.EQ.0 ) GOTO 6

C     Wire source loop.

         DO 5 J = 1 , N
            CALL TRIO(J,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            DO 3 I = 1 , JSNO
               IJ = JCO(I)
               JCO(I) = ((IJ-1)/NP)*MP2 + IJ
    3       CONTINUE
            IF ( I1.LE.IN2 ) CALL CMWW(J,I1,IN2,CCM,NROW,CCM,NROW,1,
     &                                 LD,SALP,ICON1,ICON2,X,Y,Z,SI,
     &                                 BI,CAB,SAB,JMAX,JSNO,JCO,AX,BX,
     &                                 CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMWS(J,IM1,IM2,CCM(1,IST),NROW,CCM,
     &                                  NROW,1,LD,SALP,X,Y,Z,SI,BI,CAB,
     &                                  SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &                                  JMAX,JSNO,JCO,AX,BX,CX,DEBUG)
            IF ( NLOAD.EQ.0 ) GOTO 5

C     Matrix elements modified by loading.

            IF ( J.GT.NP ) GOTO 5
            IPR = J - ISV
            IF ( IPR.LT.1 .OR. IPR.GT.IT ) GOTO 5
            ZAJ = ZARRAY(J)
            DO 4 I = 1 , JSNO
               JSS = JCO(I)
               CCM(JSS,IPR) = CCM(JSS,IPR) - (AX(I)+CX(I))*ZAJ
    4       CONTINUE
    5    CONTINUE
    6    IF ( M.EQ.0 ) GOTO 8
C     Matrix elements for patch current sources.
         JM1 = 1 - MP
         JM2 = 0
         JST = 1 - MP2
         DO 7 I = 1 , NOP
            JM1 = JM1 + MP
            JM2 = JM2 + MP
            JST = JST + NNPEQ
            IF ( I1.LE.IN2 ) CALL CMSW(JM1,JM2,I1,IN2,CCM(JST,1),CCM,0,
     &                                 NROW,1,LD,SALP,N1,N,NP,M1,M,MP,
     &                                 ICON1,ICON2,ICONX,X,Y,Z,SI,BI,
     &                                 CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &                                 JMAX,JSNO,JCO,AX,BX,CX,IFAIL,
     &                                 DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMSS(JM1,JM2,IM1,IM2,CCM(JST,IST),
     &                                  NROW,1,LD,SALP,X,Y,Z,BI,T1X,T1Y,
     &                                  T1Z,T2X,T2Y,T2Z,DEBUG)
    7    CONTINUE
    8    IF ( ICASE.EQ.1 ) GOTO 16
         IF ( ICASE.EQ.3 ) GOTO 15
C     Combine elements for symmetry modes.
         DO 14 I = 1 , IT
            DO 13 J = 1 , NNPEQ
               DO 9 K = 1 , NOP
                  KA = J + (K-1)*NNPEQ
                  D(K) = CCM(KA,I)
    9          CONTINUE
               DETER = D(1)
               DO 10 KK = 2 , NOP
                  DETER = DETER + D(KK)
   10          CONTINUE
               CCM(J,I) = DETER
               DO 12 K = 2 , NOP
                  KA = J + (K-1)*NNPEQ
                  DETER = D(1)
                  DO 11 KK = 2 , NOP
                     DETER = DETER + D(KK)*SSX(K,KK)
   11             CONTINUE
                  CCM(KA,I) = DETER
   12          CONTINUE
   13       CONTINUE
   14    CONTINUE
         IF ( ICASE.LT.3 ) GOTO 16
C     Write block for out-of-core cases.
   15    CALL BLCKOT(CCM,11,1,IOUT,DEBUG)
   16 CONTINUE
      IF ( ICASE.GT.2 ) REWIND CHTMP1
 
      RETURN
 
      END
