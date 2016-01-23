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

      SUBROUTINE CONECT( IGND , LD , SALP , IPSYM , N1 , N2 , N , NP , 
     &                   M1 , M2 , M , MP , ICON1 , ICON2 , ICONX , X , 
     &                   Y , Z , BI , X2 , Y2 , Z2 , T1X , T1Y , T1Z , 
     &                   T2X , T2Y , T2Z , JMAX , NSMAXX , NPMAX , 
     &                   NSCON , NPCON , JCO , ISCON , IPCON , IFAIL ,
     &                   DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* CONECT sets up segment connection data in arrays ICON1 and ICON2   **
C* by searching for segment ends that are in contact and segment ends **
C* which contact the centre of a surface patch.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - IGND                                                      **
C* INPUT  - LD                                                        **
C* PASSED - SALP                                                      **
C* OUTPUT - IPSYM                                                     **
C* INPUT  - N1                                                        **
C* INPUT  - N2                                                        **
C* INPUT  - N                                                         **
C* OUTPUT - NP                                                        **
C* INPUT  - M1                                                        **
C* INPUT  - M2                                                        **
C* INPUT  - M                                                         **
C* OUTPUT - MP                                                        **
C* OUTPUT - ICON1                                                     **
C* OUTPUT - ICON2                                                     **
C* OUTPUT - ICONX                                                     **
C* OUTPUT - X                                                         **
C* OUTPUT - Y                                                         **
C* OUTPUT - Z                                                         **
C* PASSED - BI                                                        **
C* OUTPUT - X2                                                        **
C* OUTPUT - Y2                                                        **
C* OUTPUT - Z2                                                        **
C* PASSED - T1X                                                       **
C* PASSED - T1Y                                                       **
C* PASSED - T1Z                                                       **
C* PASSED - T2X                                                       **
C* PASSED - T2Y                                                       **
C* PASSED - T2Z                                                       **
C* INPUT  - JMAX                                                      **
C* INPUT  - NSMAXX                                                    **
C* INPUT  - NPMAX                                                     **
C* OUTPUT - NSCON                                                     **
C* OUTPUT - NPCON                                                     **
C* OUTPUT - JCO                                                       **
C* OUTPUT - ISCON                                                     **
C* OUTPUT - IPCON                                                     **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       SUBPH                                                  **
C* called by   DATAGN                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL SUBPH

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER IGND , LD , IPSYM , N1 , N2 , N , NP , M1 , M2 , M , MP , 
     &        JMAX , NSMAXX , NPMAX , NSCON , NPCON , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , ICONX(LD) , 
     &        JCO(JMAX) , ISCON(NSMAXX) , IPCON(NPMAX)
      REAL*8 SALP(LD) , X(LD) , Y(LD) , Z(LD) , 
     &       BI(LD) , X2(LD) , Y2(LD) , Z2(LD) , 
     &       T1X(LD) , T1Y(LD) , T1Z(LD) , T2X(LD) , 
     &       T2Y(LD) , T2Z(LD)

C     Local variables.
      INTEGER I , IC , IEND , ISEG , IX , J , JEND , NSFLG
      REAL*8 SEP , SLEN , SMIN , XA , XI1 , XI2 , XXS , YA, YI1 , YI2 , 
     &       YYS , ZA , ZI1 , ZI2 , ZZS
 
C     Data initialisation.
      DATA SMIN/1.0D-3/
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING CONECT'
      ENDIF
 
      NSCON = 0
      NPCON = 0
      IF ( IGND.EQ.0 ) GOTO 3
      WRITE (CHRSLT,54)
      IF ( IGND.GT.0 ) WRITE (CHRSLT,55)
      IF ( IPSYM.NE.2 ) GOTO 1
      NP = 2*NP
      MP = 2*MP
    1 IF ( IABS(IPSYM).LE.2 ) GOTO 2
      NP = N
      MP = M
    2 IF ( NP.GT.N ) THEN
         IFAIL=4
         RETURN
      ENDIF
      IF ( NP.EQ.N .AND. MP.EQ.M ) IPSYM = 0
    3 IF ( N.EQ.0 ) GOTO 26
      DO 15 I = 1 , N
         ICONX(I) = 0
         XI1 = X(I)
         YI1 = Y(I)
         ZI1 = Z(I)
         XI2 = X2(I)
         YI2 = Y2(I)
         ZI2 = Z2(I)
         SLEN = DSQRT((XI2-XI1)**2+(YI2-YI1)**2+(ZI2-ZI1)**2)*SMIN

C     Determine connection data for end 1 of segment.

         IF ( IGND.LT.1 ) GOTO 5
         IF ( ZI1.GT.-SLEN ) GOTO 4
         WRITE (CHRSLT,56) I
         IFAIL=5
         RETURN
    4    IF ( ZI1.GT.SLEN ) GOTO 5
         ICON1(I) = I
         Z(I) = 0.0D0
         GOTO 9
    5    IC = I
         DO 7 J = 2 , N
            IC = IC + 1
            IF ( IC.GT.N ) IC = 1
            SEP = DABS(XI1-X(IC)) + DABS(YI1-Y(IC)) + DABS(ZI1-Z(IC))
            IF ( SEP.GT.SLEN ) GOTO 6
            ICON1(I) = -IC
            GOTO 8
    6       SEP = DABS(XI1-X2(IC)) + DABS(YI1-Y2(IC)) + DABS(ZI1-Z2(IC))
            IF ( SEP.GT.SLEN ) GOTO 7
            ICON1(I) = IC
            GOTO 8
    7    CONTINUE
         IF ( I.LT.N2 .AND. ICON1(I).GT.PATOFF ) GOTO 8
         ICON1(I) = 0

C     Determine connection data for end 2 of segment.

    8    IF ( IGND.LT.1 ) GOTO 12
    9    IF ( ZI2.GT.-SLEN ) GOTO 10
         WRITE (CHRSLT,56) I
         IFAIL=5
         RETURN
   10    IF ( ZI2.GT.SLEN ) GOTO 12
         IF ( ICON1(I).NE.I ) GOTO 11
         WRITE (CHRSLT,57) I
         IFAIL=6
         RETURN
   11    ICON2(I) = I
         Z2(I) = 0.0D0
         GOTO 15
   12    IC = I
         DO 14 J = 2 , N
            IC = IC + 1
            IF ( IC.GT.N ) IC = 1
            SEP = DABS(XI2-X(IC)) + DABS(YI2-Y(IC)) + DABS(ZI2-Z(IC))
            IF ( SEP.GT.SLEN ) GOTO 13
            ICON2(I) = IC
            GOTO 15
   13       SEP = DABS(XI2-X2(IC)) + DABS(YI2-Y2(IC)) + DABS(ZI2-Z2(IC))
            IF ( SEP.GT.SLEN ) GOTO 14
            ICON2(I) = -IC
            GOTO 15
   14    CONTINUE
         IF ( I.LT.N2 .AND. ICON2(I).GT.PATOFF ) GOTO 15
         ICON2(I) = 0
   15 CONTINUE
      IF ( M.EQ.0 ) GOTO 26
C     Find wire-surface connections for new patches.
      IX = LD + 1 - M1
      I = M2
   16 IF ( I.GT.M ) GOTO 20
      IX = IX - 1
      XXS = X(IX)
      YYS = Y(IX)
      ZZS = Z(IX)
      DO 18 ISEG = 1 , N
         XI1 = X(ISEG)
         YI1 = Y(ISEG)
         ZI1 = Z(ISEG)
         XI2 = X2(ISEG)
         YI2 = Y2(ISEG)
         ZI2 = Z2(ISEG)
         SLEN = (DABS(XI2-XI1)+DABS(YI2-YI1)+DABS(ZI2-ZI1))*SMIN
C     For first end of segment.
         SEP = DABS(XI1-XXS) + DABS(YI1-YYS) + DABS(ZI1-ZZS)
         IF ( SEP.GT.SLEN ) GOTO 17
C     Connection - divide patch into 4 patches at present array loc.
         ICON1(ISEG) = PATOFF + I
         IC = 0
         CALL SUBPH(I,IC,LD,SALP,M,MP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &              DEBUG)
         GOTO 19
   17    SEP = DABS(XI2-XXS) + DABS(YI2-YYS) + DABS(ZI2-ZZS)
         IF ( SEP.GT.SLEN ) GOTO 18
         ICON2(ISEG) = PATOFF + I
         IC = 0
         CALL SUBPH(I,IC,LD,SALP,M,MP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &              DEBUG)
         GOTO 19
   18 CONTINUE
   19 I = I + 1
      GOTO 16
C     Repeat search for new segments connected to ngf patches.
   20 IF ( M1.EQ.0 .OR. N2.GT.N ) GOTO 26
      IX = LD + 1
      I = 1
   21 IF ( I.GT.M1 ) GOTO 25
      IX = IX - 1
      XXS = X(IX)
      YYS = Y(IX)
      ZZS = Z(IX)
      DO 23 ISEG = N2 , N
         XI1 = X(ISEG)
         YI1 = Y(ISEG)
         ZI1 = Z(ISEG)
         XI2 = X2(ISEG)
         YI2 = Y2(ISEG)
         ZI2 = Z2(ISEG)
         SLEN = (DABS(XI2-XI1)+DABS(YI2-YI1)+DABS(ZI2-ZI1))*SMIN
         SEP = DABS(XI1-XXS) + DABS(YI1-YYS) + DABS(ZI1-ZZS)
         IF ( SEP.GT.SLEN ) GOTO 22
         ICON1(ISEG) = PATOFF + 1 + M
         IC = 1
         NPCON = NPCON + 1
         IPCON(NPCON) = I
         CALL SUBPH(I,IC,LD,SALP,M,MP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &              DEBUG)
         GOTO 24
   22    SEP = DABS(XI2-XXS) + DABS(YI2-YYS) + DABS(ZI2-ZZS)
         IF ( SEP.GT.SLEN ) GOTO 23
         ICON2(ISEG) = PATOFF + 1 + M
         IC = 1
         NPCON = NPCON + 1
         IPCON(NPCON) = I
         CALL SUBPH(I,IC,LD,SALP,M,MP,X,Y,Z,BI,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
     &              DEBUG)
         GOTO 24
   23 CONTINUE
   24 I = I + 1
      GOTO 21
   25 IF ( NPCON.LE.NPMAX ) GOTO 26
      WRITE (CHRSLT,62) NPMAX
      IFAIL=7
      RETURN
   26 WRITE (CHRSLT,58) N , NP , IPSYM
      IF ( M.GT.0 ) WRITE (CHRSLT,61) M , MP
      ISEG = (N+M)/(NP+MP)
      IF ( ISEG.EQ.1 ) GOTO 30
      IF ( IPSYM ) 28 , 27 , 29
   27 IFAIL=8
      RETURN
   28 WRITE (CHRSLT,59) ISEG
      GOTO 30
   29 IC = ISEG/2
      IF ( ISEG.EQ.8 ) IC = 3
      WRITE (CHRSLT,60) IC
   30 IF ( N.EQ.0 ) GOTO 48
      WRITE (CHRSLT,50)
      ISEG = 0
C     Adjust connected seg. ends to exactly coincide. Print junctions
C     of 3 or more seg. Also find old seg. connecting to new seg.
      DO 44 J = 1 , N
         IEND = -1
         JEND = -1
         IX = ICON1(J)
         IC = 1
         JCO(1) = -J
         XA = X(J)
         YA = Y(J)
         ZA = Z(J)
   31    IF ( IX.EQ.0 ) GOTO 43
         IF ( IX.EQ.J ) GOTO 43
         IF ( IX.GT.PATOFF ) GOTO 43
         NSFLG = 0
   32    IF ( IX ) 33 , 49 , 34
   33    IX = -IX
         GOTO 35
   34    JEND = -JEND
   35    IF ( IX.EQ.J ) GOTO 37
         IF ( IX.LT.J ) GOTO 43
         IC = IC + 1
         IF ( IC.GT.JMAX ) GOTO 49
         JCO(IC) = IX*JEND
         IF ( IX.GT.N1 ) NSFLG = 1
         IF ( JEND.EQ.1 ) GOTO 36
         XA = XA + X(IX)
         YA = YA + Y(IX)
         ZA = ZA + Z(IX)
         IX = ICON1(IX)
         GOTO 32
   36    XA = XA + X2(IX)
         YA = YA + Y2(IX)
         ZA = ZA + Z2(IX)
         IX = ICON2(IX)
         GOTO 32
   37    SEP = IC
         XA = XA/SEP
         YA = YA/SEP
         ZA = ZA/SEP
         DO 39 I = 1 , IC
            IX = JCO(I)
            IF ( IX.GT.0 ) GOTO 38
            IX = -IX
            X(IX) = XA
            Y(IX) = YA
            Z(IX) = ZA
            GOTO 39
   38       X2(IX) = XA
            Y2(IX) = YA
            Z2(IX) = ZA
   39    CONTINUE
         IF ( N1.EQ.0 ) GOTO 42
         IF ( NSFLG.EQ.0 ) GOTO 42
         DO 41 I = 1 , IC
            IX = IABS(JCO(I))
            IF ( IX.GT.N1 ) GOTO 41
            IF ( ICONX(IX).NE.0 ) GOTO 41
            NSCON = NSCON + 1
            IF ( NSCON.LE.NSMAXX ) GOTO 40
            WRITE (CHRSLT,62) NSMAXX
            IFAIL=7
            RETURN
   40       ISCON(NSCON) = IX
            ICONX(IX) = NSCON
   41    CONTINUE
   42    IF ( IC.LT.3 ) GOTO 43
         ISEG = ISEG + 1
         WRITE (CHRSLT,51) ISEG , (JCO(I),I=1,IC)
   43    IF ( IEND.EQ.1 ) GOTO 44
         IEND = 1
         JEND = 1
         IX = ICON2(J)
         IC = 1
         JCO(1) = J
         XA = X2(J)
         YA = Y2(J)
         ZA = Z2(J)
         GOTO 31
   44 CONTINUE
      IF ( ISEG.EQ.0 ) WRITE (CHRSLT,52)
      IF ( N1.EQ.0 .OR. M1.EQ.M ) GOTO 48
C     Find old segments that connect to new patches.
      DO 47 J = 1 , N1
         IX = ICON1(J)
         IF ( IX.LT.PATOFF ) GOTO 45
         IX = IX - PATOFF
         IF ( IX.GT.M1 ) GOTO 46
   45    IX = ICON2(J)
         IF ( IX.LT.PATOFF ) GOTO 47
         IX = IX - PATOFF
         IF ( IX.LT.M2 ) GOTO 47
   46    IF ( ICONX(J).NE.0 ) GOTO 47
         NSCON = NSCON + 1
         ISCON(NSCON) = J
         ICONX(J) = NSCON
   47 CONTINUE
   48 CONTINUE
 
      RETURN
 
   49 WRITE (CHRSLT,53) IX
 
      IFAIL=9
      RETURN
C
   50 FORMAT (//,9X,'- MULTIPLE WIRE JUNCTIONS -',/,1X,'JUNCTION',4X,
     &        'SEGMENTS  (- FOR END 1, + FOR END 2)')
   51 FORMAT (1X,I5,5X,20I5,/,(11X,20I5))
   52 FORMAT (2X,'NONE')
   53 FORMAT (' CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT',I5)
   54 FORMAT (/,3X,'GROUND PLANE SPECIFIED.')
   55 FORMAT (/,3X,'WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL BE ',
     &        'INTERPOLATED TO IMAGE IN GROUND PLANE.',/)
   56 FORMAT (' GEOMETRY DATA ERROR-- SEGMENT',I5,
     &        ' EXTENDS BELOW GROUND')
   57 FORMAT (' GEOMETRY DATA ERROR--SEGMENT',I5,' LIES IN GROUND ',
     &        'PLANE.')
   58 FORMAT (/,3X,'TOTAL SEGMENTS USED=',I5,5X,'NO. SEG. IN ',
     &        'A SYMMETRIC CELL=',I5,5X,'SYMMETRY FLAG=',I3)
   59 FORMAT (' STRUCTURE HAS',I4,' FOLD ROTATIONAL SYMMETRY',/)
   60 FORMAT (' STRUCTURE HAS',I2,' PLANES OF SYMMETRY',/)
   61 FORMAT (3X,'TOTAL PATCHES USED=',I5,6X,
     &        'NO. PATCHES IN A SYMMETRIC CELL=',I5)
   62 FORMAT (
     &' ERROR - NO. NEW SEGMENTS CONNECTED TO N.G.F. SEGMENTSOR PATCHES 
     &EXCEEDS LIMIT OF,',I5)
 
      END
