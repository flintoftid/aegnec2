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

      SUBROUTINE CMNGF( CB , CC , CD , NB , NC , ND , RKHX , IEXKX , 
     &                  LD , SALP , NLOAD , NLODF , ZARRAY , N1 , N2 , 
     &                  N , NP , M1 , M2 , M , MP , ICON1 , ICON2 , 
     &                  ICONX , X , Y , Z , SI , BI , CAB , SAB , T1X , 
     &                  T1Y , T1Z , T2X , T2Y , T2Z , JMAX , NSMAXX , 
     &                  NPMAX , JSNO , NSCON , NPCON , JCO , ISCON , 
     &                  IPCON , AX , BX , CX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CMNGF fills interaction matricies B, C, and D for NGF solution.    **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - CB                                                        **
C* OUTPUT - CC                                                        **
C* OUTPUT - CD                                                        **
C* INPUT  - NB                                                        **
C* INPUT  - NC                                                        **
C* INPUT  - ND                                                        **
C* INPUT  - RKHX                                                      **
C* INPUT  - IEXKX                                                     **
C* PASSED - LD                                                        **
C* PASSED - SALP                                                      **
C* INPUT  - NLOAD                                                     **
C* INPUT  - NLODF                                                     **
C* INPUT  - ZARRAY                                                    **
C* INPUT  - N1                                                        **
C* INPUT  - N2                                                        **
C* INPUT  - N                                                         **
C* PASSED - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M2                                                        **
C* INPUT  - M                                                         **
C* PASSED - MP                                                        **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - ICONX                                                     **
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
C* INPUT  - NSMAXX                                                    **
C* INPUT  - NPMAX                                                     **
C* OUTPUT - JSNO                                                      **
C* INPUT  - NSCON                                                     **
C* INPUT  - NPCON                                                     **
C* OUTPUT - JCO                                                       **
C* INPUT  - ISCON                                                     **
C* INPUT  - IPCON                                                     **
C* OUTPUT - AX                                                        **
C* OUTPUT - BX                                                        **
C* OUTPUT - CX                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXK     IEXK    RKH                                    **
C* uses value  EXK     ICASX   NBBL    NBBX    NLBL    NLBX    NPBL   **
C*             NPBX                                                   **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       CMSS    CMSW    CMWS    CMWW    TBF     TRIO           **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL CMSS , CMSW , CMWS , CMWW , TBF , TRIO

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER NB , NC , ND , IEXKX , LD , NLOAD , NLODF , N1 , N2 , N , 
     &        NP , M1 , M2 , M , MP , JMAX , NSMAXX , NPMAX , JSNO , 
     &        NSCON , NPCON , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , ICONX(LD) , 
     &        JCO(JMAX) , ISCON(NSMAXX) , IPCON(NPMAX)
      REAL*8 RKHX
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       SI(LD) ,  BI(LD) ,CAB(LD) , SAB(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD) ,AX(JMAX) , BX(JMAX) , CX(JMAX) 
      COMPLEX*16 CB(NB,1) , CC(NC,1) , CD(ND,1) , ZARRAY(LD)

C     Local variables.
      INTEGER I , I1 , I2 , IBLK , IM1 , IM2 , IMX , IN1 , IN2 ,
     &        IR , IST , ISV , ISVV , IT , ITX , IX , J , JSS , JSX , 
     &        JX , M1EQ , M2EQ , MEQ , NEQN , NEQP , NEQS , NEQSP

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'matpar.inc'

      
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CMNGF'
      ENDIF
      
      RKH = RKHX
      IEXK = IEXKX
      M1EQ = 2*M1
      M2EQ = M1EQ + 1
      MEQ = 2*M
      NEQP = ND - NPCON*2
      NEQS = NEQP - NSCON
      NEQSP = NEQS + NC
      NEQN = NC + N - N1
      ITX = 1
      IF ( NSCON.GT.0 ) ITX = 2
      IF ( ICASX.EQ.1 ) GOTO 1
      REWIND CHTMP2
      REWIND CHTMP4
      REWIND CHTMP5
      IF ( ICASX.GT.2 ) GOTO 5
    1 DO 4 J = 1 , ND
         DO 2 I = 1 , ND
            CD(I,J) = (0.0D0,0.0D0)
    2    CONTINUE
         DO 3 I = 1 , NB
            CB(I,J) = (0.0D0,0.0D0)
            CC(I,J) = (0.0D0,0.0D0)
    3    CONTINUE
    4 CONTINUE
    5 IST = N - N1 + 1
      IT = NPBX
      ISV = -NPBX
C     Loop thru 24 fills B. For ICASX=1 or 2 also fills D(WW), D(WS).
      DO 26 IBLK = 1 , NBBX
         ISV = ISV + NPBX
         IF ( IBLK.EQ.NBBX ) IT = NLBX
         IF ( ICASX.LT.3 ) GOTO 8
         DO 7 J = 1 , ND
            DO 6 I = 1 , IT
               CB(I,J) = (0.0D0,0.0D0)
    6       CONTINUE
    7    CONTINUE
    8    I1 = ISV + 1
         I2 = ISV + IT
         IN2 = I2
         IF ( IN2.GT.N1 ) IN2 = N1
         IM1 = I1 - N1
         IM2 = I2 - N1
         IF ( IM1.LT.1 ) IM1 = 1
         IMX = 1
         IF ( I1.LE.N1 ) IMX = N1 - I1 + 2
         IF ( N2.GT.N ) GOTO 13
C     Fill B(WW),B(WS). For ICASX=1,2 fill D(WW),D(WS).
         DO 12 J = N2 , N
            CALL TRIO(J,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            DO 10 I = 1 , JSNO
               JSS = JCO(I)
               IF ( JSS.LT.N2 ) GOTO 9
C     Set JCO when source is new basis function on new segment.
               JCO(I) = JSS - N1
               GOTO 10
C     Source is portion of modified basis function on new segment.
    9          JCO(I) = NEQS + ICONX(JSS)
   10       CONTINUE
            IF ( I1.LE.IN2 ) CALL CMWW(J,I1,IN2,CB,NB,CB,NB,0,LD,
     &                                 SALP,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &                                 SAB,JMAX,JSNO,JCO,AX,BX,CX,IFAIL,
     &                                 DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMWS(J,IM1,IM2,CB(IMX,1),NB,CB,NB,0,
     &                                  LD,SALP,X,Y,Z,SI,BI,CAB,SAB,T1X,
     &                                  T1Y,T1Z,T2X,T2Y,T2Z,JMAX,JSNO,
     &                                  JCO,AX,BX,CX,DEBUG)
            IF ( ICASX.GT.2 ) GOTO 12
            CALL CMWW(J,N2,N,CD,ND,CD,ND,1,LD,SALP,ICON1,ICON2,X,Y,
     &                Z,SI,BI,CAB,SAB,JMAX,JSNO,JCO,AX,BX,CX,IFAIL,
     &                DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( M2.LE.M ) CALL CMWS(J,M2EQ,MEQ,CD(1,IST),ND,CD,ND,1,
     &                               LD,SALP,X,Y,Z,SI,BI,CAB,SAB,T1X,
     &                               T1Y,T1Z,T2X,T2Y,T2Z,JMAX,JSNO,JCO,
     &                               AX,BX,CX,DEBUG)
C     Loading in D(WW).
            IF ( NLOAD.EQ.0 ) GOTO 12
            IR = J - N1
            EXK = ZARRAY(J)
            DO 11 I = 1 , JSNO
               JSS = JCO(I)
               CD(JSS,IR) = CD(JSS,IR) - (AX(I)+CX(I))*EXK
   11       CONTINUE
   12    CONTINUE
   13    IF ( NSCON.EQ.0 ) GOTO 22
C     Fill  B(WW) prime.
         DO 21 I = 1 , NSCON
            J = ISCON(I)
C     Sources are new or modified basis functions on old segments which
C     connect to new segments.
            CALL TRIO(J,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            JSS = 0
            DO 16 IX = 1 , JSNO
               IR = JCO(IX)
               IF ( IR.LT.N2 ) GOTO 14
               IR = IR - N1
               GOTO 15
   14          IR = ICONX(IR)
               IF ( IR.EQ.0 ) GOTO 16
               IR = NEQS + IR
   15          JSS = JSS + 1
               JCO(JSS) = IR
               AX(JSS) = AX(IX)
               BX(JSS) = BX(IX)
               CX(JSS) = CX(IX)
   16       CONTINUE
            JSNO = JSS
            IF ( I1.LE.IN2 ) CALL CMWW(J,I1,IN2,CB,NB,CB,NB,0,LD,
     &                                 SALP,ICON1,ICON2,X,Y,Z,SI,BI,CAB,
     &                                 SAB,JMAX,JSNO,JCO,AX,BX,CX,IFAIL,
     &                                 DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMWS(J,IM1,IM2,CB(IMX,1),NB,CB,NB,0,
     &                                  LD,SALP,X,Y,Z,SI,BI,CAB,SAB,
     &                                  T1X,T1Y,T1Z,T2X,T2Y,T2Z,JMAX,
     &                                  JSNO,JCO,AX,BX,CX,DEBUG)
C     Source is singular component of patch current that is part of
C     modified basis function for old segment that connects to a new
C     segment on end opposite patch.
            IF ( I1.LE.IN2 ) CALL CMSW(J,I,I1,IN2,CB,CB,0,NB,-1,
     &                                 LD,SALP,N1,N,NP,M1,M,MP,ICON1,
     &                                 ICON2,ICONX,X,Y,Z,SI,BI,CAB,SAB,
     &                                 T1X,T1Y,T1Z,T2X,T2Y,T2Z,JMAX,
     &                                 JSNO,JCO,AX,BX,CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( NLODF.EQ.0 ) GOTO 18
            JX = J - ISV
            IF ( JX.LT.1 .OR. JX.GT.IT ) GOTO 18
            EXK = ZARRAY(J)
            DO 17 IX = 1 , JSNO
               JSS = JCO(IX)
               CB(JX,JSS) = CB(JX,JSS) - (AX(IX)+CX(IX))*EXK
   17       CONTINUE
C     Sources are portions of modified basis function J on old segments
C     excluding old segments that directly connect to new segments.
   18       CALL TBF(J,1,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &               CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            JSX = JSNO
            JSNO = 1
            IR = JCO(1)
            JCO(1) = NEQS + I
            DO 20 IX = 1 , JSX
               IF ( IX.EQ.1 ) GOTO 19
               IR = JCO(IX)
               AX(1) = AX(IX)
               BX(1) = BX(IX)
               CX(1) = CX(IX)
   19          IF ( IR.GT.N1 ) GOTO 20
               IF ( ICONX(IR).NE.0 ) GOTO 20
               IF ( I1.LE.IN2 ) CALL CMWW(IR,I1,IN2,CB,NB,CB,NB,0,
     &                                    LD,SALP,ICON1,ICON2,X,Y,Z,
     &                                    SI,BI,CAB,SAB,JMAX,JSNO,JCO,
     &                                    AX,BX,CX,IFAIL,DEBUG)
               IF(IFAIL.NE.0) RETURN
               IF ( IM1.LE.IM2 ) CALL CMWS(IR,IM1,IM2,CB(IMX,1),NB,CB,
     &                                     NB,0,LD,SALP,X,Y,Z,SI,BI,CAB,
     &                                     SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &                                     JMAX,JSNO,JCO,AX,BX,CX,DEBUG)
C     Loading for B(WW)prime.
               IF ( NLODF.EQ.0 ) GOTO 20
               JX = IR - ISV
               IF ( JX.LT.1 .OR. JX.GT.IT ) GOTO 20
               EXK = ZARRAY(IR)
               JSS = JCO(1)
               CB(JX,JSS) = CB(JX,JSS) - (AX(1)+CX(1))*EXK
   20       CONTINUE
   21    CONTINUE
   22    IF ( NPCON.EQ.0 ) GOTO 24
         JSS = NEQP
C     Fill B(SS)prime to set old patch basis functions to zero for
c     patches that connect to new segments.
         DO 23 I = 1 , NPCON
            IX = IPCON(I)*2 + N1 - ISV
            IR = IX - 1
            JSS = JSS + 1
            IF ( IR.GT.0 .AND. IR.LE.IT ) CB(IR,JSS) = (1.0D0,0.0D0)
            JSS = JSS + 1
            IF ( IX.GT.0 .AND. IX.LE.IT ) CB(IX,JSS) = (1.0D0,0.0D0)
   23    CONTINUE
   24    IF ( M2.GT.M ) GOTO 25
C     Fill B(SW) and B(SS)
         IF ( I1.LE.IN2 ) CALL CMSW(M2,M,I1,IN2,CB(1,IST),CB,N1,NB,0,
     &                              LD,SALP,N1,N,NP,M1,M,MP,ICON1,ICON2,
     &                              ICONX,X,Y,Z,SI,BI,CAB,SAB,T1X,T1Y,
     &                              T1Z,T2X,T2Y,T2Z,JMAX,JSNO,JCO,AX,BX,
     &                              CX,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF ( IM1.LE.IM2 ) CALL CMSS(M2,M,IM1,IM2,CB(IMX,IST),NB,0,
     &                               LD,SALP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,
     &                               T2Y,T2Z,DEBUG)
   25    IF ( ICASX.EQ.1 ) GOTO 26
         WRITE (CHTMP4) ((CB(I,J),I=1,IT),J=1,ND)
   26 CONTINUE
C     Filling  B complete. Start on C and D.
      IT = NPBL
      ISV = -NPBL
      DO 46 IBLK = 1 , NBBL
         ISV = ISV + NPBL
         ISVV = ISV + NC
         IF ( IBLK.EQ.NBBL ) IT = NLBL
         IF ( ICASX.LT.3 ) GOTO 30
         DO 29 J = 1 , IT
            DO 27 I = 1 , NC
               CC(I,J) = (0.0D0,0.0D0)
   27       CONTINUE
            DO 28 I = 1 , ND
               CD(I,J) = (0.0D0,0.0D0)
   28       CONTINUE
   29    CONTINUE
   30    I1 = ISVV + 1
         I2 = ISVV + IT
         IN1 = I1 - M1EQ
         IN2 = I2 - M1EQ
         IF ( IN2.GT.N ) IN2 = N
         IM1 = I1 - N
         IM2 = I2 - N
         IF ( IM1.LT.M2EQ ) IM1 = M2EQ
         IF ( IM2.GT.MEQ ) IM2 = MEQ
         IMX = 1
         IF ( IN1.LE.IN2 ) IMX = NEQN - I1 + 2
         IF ( ICASX.LT.3 ) GOTO 35
         IF ( N2.GT.N ) GOTO 35
C     Same as DO 24 loop to fill D(WW) for ICASX greater than 2.
         DO 34 J = N2 , N
            CALL TRIO(J,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            DO 32 I = 1 , JSNO
               JSS = JCO(I)
               IF ( JSS.LT.N2 ) GOTO 31
               JCO(I) = JSS - N1
               GOTO 32
   31          JCO(I) = NEQS + ICONX(JSS)
   32       CONTINUE
            IF ( IN1.LE.IN2 ) CALL CMWW(J,IN1,IN2,CD,ND,CD,ND,1,LD,
     &                                  SALP,ICON1,ICON2,X,Y,Z,SI,BI,
     &                                  CAB,SAB,JMAX,JSNO,JCO,AX,BX,CX,
     &                                  IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMWS(J,IM1,IM2,CD(1,IMX),ND,CD,ND,1,
     &                                  LD,SALP,X,Y,Z,SI,BI,CAB,SAB,T1X,
     &                                  T1Y,T1Z,T2X,T2Y,T2Z,JMAX,JSNO,
     &                                  JCO,AX,BX,CX,DEBUG)
            IF ( NLOAD.EQ.0 ) GOTO 34
            IR = J - N1 - ISV
            IF ( IR.LT.1 .OR. IR.GT.IT ) GOTO 34
            EXK = ZARRAY(J)
            DO 33 I = 1 , JSNO
               JSS = JCO(I)
               CD(JSS,IR) = CD(JSS,IR) - (AX(I)+CX(I))*EXK
   33       CONTINUE
   34    CONTINUE
   35    IF ( M2.GT.M ) GOTO 36
C     Fill D(SW) and D(SS).
         IF ( IN1.LE.IN2 ) CALL CMSW(M2,M,IN1,IN2,CD(IST,1),CD,N1,ND,1,
     &                               LD,SALP,N1,N,NP,M1,M,MP,ICON1,
     &                               ICON2,ICONX,X,Y,Z,SI,BI,CAB,SAB,
     &                               T1X,T1Y,T1Z,T2X,T2Y,T2Z,JMAX,JSNO,
     &                               JCO,AX,BX,CX,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF ( IM1.LE.IM2 ) CALL CMSS(M2,M,IM1,IM2,CD(IST,IMX),ND,1,
     &                               LD,SALP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,
     &                               T2Y,T2Z,DEBUG)
   36    IF ( N1.LT.1 ) GOTO 42
C     Fill C(WW),C(WS), D(WW)prime, and D(WS)prime.
         DO 40 J = 1 , N1
            CALL TRIO(J,LD,ICON1,ICON2,SI,BI,JMAX,JSNO,JCO,AX,BX,
     &                CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( NSCON.EQ.0 ) GOTO 39
            DO 38 IX = 1 , JSNO
               JSS = JCO(IX)
               IF ( JSS.LT.N2 ) GOTO 37
               JCO(IX) = JSS + M1EQ
               GOTO 38
   37          IR = ICONX(JSS)
               IF ( IR.NE.0 ) JCO(IX) = NEQSP + IR
   38       CONTINUE
   39       IF ( IN1.LE.IN2 ) CALL CMWW(J,IN1,IN2,CC,NC,CD,ND,ITX,
     &                                  LD,SALP,ICON1,ICON2,X,Y,Z,
     &                                  SI,BI,CAB,SAB,JMAX,JSNO,JCO,AX,
     &                                  BX,CX,IFAIL,DEBUG)
            IF(IFAIL.NE.0) RETURN
            IF ( IM1.LE.IM2 ) CALL CMWS(J,IM1,IM2,CC(1,IMX),NC,
     &                                  CD(1,IMX),ND,ITX,LD,SALP,X,Y,Z,
     &                                  SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,
     &                                  T2Y,T2Z,JMAX,JSNO,JCO,AX,BX,CX,
     &                                  DEBUG)
   40    CONTINUE
         IF ( NSCON.EQ.0 ) GOTO 42
C     Fill C(WW)prime.
         DO 41 IX = 1 , NSCON
            IR = ISCON(IX)
            JSS = NEQS + IX - ISV
            IF ( JSS.GT.0 .AND. JSS.LE.IT ) CC(IR,JSS) = (1.0D0,0.0D0)
   41    CONTINUE
   42    IF ( NPCON.EQ.0 ) GOTO 44
         JSS = NEQP - ISV
C     Fill C(SS)prime.
         DO 43 I = 1 , NPCON
            IX = IPCON(I)*2 + N1
            IR = IX - 1
            JSS = JSS + 1
            IF ( JSS.GT.0 .AND. JSS.LE.IT ) CC(IR,JSS) = (1.0D0,0.0D0)
            JSS = JSS + 1
            IF ( JSS.GT.0 .AND. JSS.LE.IT ) CC(IX,JSS) = (1.0D0,0.0D0)
   43    CONTINUE
   44    IF ( M1.LT.1 ) GOTO 45
C     Fill C(SW) and C(SS).
         IF ( IN1.LE.IN2 ) CALL CMSW(1,M1,IN1,IN2,CC(N2,1),CC,0,NC,1,
     &                               LD,SALP,N1,N,NP,M1,M,MP,ICON1,
     &                               ICON2,ICONX,X,Y,Z,SI,BI,CAB,SAB,
     &                               T1X,T1Y,T1Z,T2X,T2Y,T2Z,JMAX,JSNO,
     &                               JCO,AX,BX,CX,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF ( IM1.LE.IM2 ) CALL CMSS(1,M1,IM1,IM2,CC(N2,IMX),NC,1,
     &                               LD,SALP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,
     &                               T2Y,T2Z,DEBUG)
   45    CONTINUE
         IF ( ICASX.EQ.1 ) GOTO 46
         WRITE (CHTMP2) ((CD(J,I),J=1,ND),I=1,IT)
         WRITE (CHTMP5) ((CC(J,I),J=1,NC),I=1,IT)
   46 CONTINUE
      IF ( ICASX.EQ.1 ) RETURN
      REWIND CHTMP2
      REWIND CHTMP4
      REWIND CHTMP5
 
      RETURN
 
      END
