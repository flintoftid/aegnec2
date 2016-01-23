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

      SUBROUTINE NFPAT( LD , SALP , AIR , AII , BIR , BII , CIR , 
     &                  CII , CUR , N , M , ICON1 , ICON2 , WLAM , X , 
     &                  Y , Z , SI , BI , CAB , SAB , T1X , T1Y , T1Z , 
     &                  T2X , T2Y , T2Z , IPLP1 , IPLP2 , IPLP3 , 
     &                  IPLP4 , NEAR , NFEH , NRX , NRY , NRZ , XNR , 
     &                  YNR , ZNR , DXNR , DYNR , DZNR , JOBNAM , 
     &                  FILBUF , MULTIF , NRECNT , NRHCNT , IFAIL ,
     &                  DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute near E or H fields over a range of points.                 **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - LD                                                        **
C* PASSED - SALP                                                      **
C* PASSED - AIR                                                       **
C* PASSED - AII                                                       **
C* PASSED - BIR                                                       **
C* PASSED - BII                                                       **
C* PASSED - CIR                                                       **
C* PASSED - CII                                                       **
C* PASSED - CUR                                                       **
C* PASSED - N                                                         **
C* PASSED - M                                                         **
C* PASSED - ICON1                                                     **
C* PASSED - ICON2                                                     **
C* INPUT  - WLAM                                                      **
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
C* INPUT  - IPLP1                                                     **
C* INPUT  - IPLP2                                                     **
C* INPUT  - IPLP3                                                     **
C* INPUT  - IPLP4                                                     **
C* INPUT  - NEAR                                                      **
C* INPUT  - NFEH                                                      **
C* INPUT  - NRX                                                       **
C* INPUT  - NRY                                                       **
C* INPUT  - NRZ                                                       **
C* INPUT  - XNR                                                       **
C* INPUT  - YNR                                                       **
C* INPUT  - ZNR                                                       **
C* INPUT  - DXNR                                                      **
C* INPUT  - DYNR                                                      **
C* INPUT  - DZNR                                                      **
C* PASSED - JOBNAM                                                    **
C* INPUT  - FILBUF                                                    **
C* INPUT  - MULTIF                                                    **
C* OUTPUT - NRECNT                                                    **
C* OUTPUT - NRHCNT                                                    **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       CANG    FILSUF  NEFLD   NHFLD                          **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      REAL*8 CANG
      EXTERNAL CANG , NEFLD , NHFLD , FILSUF

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF
      INTEGER LD , N , M , IPLP1 , IPLP2 , IPLP3 , IPLP4 , 
     &        NEAR , NFEH , NRX , NRY , NRZ , MULTIF , NRECNT , NRHCNT , 
     &        IFAIL , DEBUG 
      INTEGER ICON1(2*LD) , ICON2(2*LD)
      REAL*8 WLAM , XNR , YNR , ZNR , DXNR , DYNR , DZNR
      REAL*8 SALP(LD) , AIR(LD) , AII(LD) , BIR(LD) , BII(LD) , 
     &       CIR(LD) , CII(LD) , X(LD) , Y(LD) , Z(LD) , SI(LD) , 
     &       BI(LD) , CAB(LD) , SAB(LD) , T1X(LD) , T1Y(LD) , 
     &       T1Z(LD) , T2X(LD) , T2Y(LD) , T2Z(LD)
      COMPLEX*16 CUR(3*LD)

C     Local variables.
      INTEGER I , J , KK
      REAL*8 AT , BT , CPH , CT , CTH , ETOTAL , EX2 , EY2 , EZ2 , 
     &       RPD , SPH , STH , TA , TMP1 , TMP2 , TMP3 , TMP4 , 
     &       TMP5 , TMP6 , XNRT , XOB , XXX , YNRT , YOB , ZNRT , ZOB
      COMPLEX*16 EX , EY , EZ
 
C     Data initialisation.
      DATA TA /0.0174532925199432957692369D0/
      DATA RPD /57.29577951308232087679815D0/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING NFPAT'
      ENDIF 
 
      IF ( NFEH.EQ.1 ) GOTO 1
      WRITE (CHRSLT,12)
      IF ( MULTIF.EQ.1 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,NRECNT,SUFNRE,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         OPEN(UNIT=CHNREF,STATUS='NEW',FILE=FILBUF)
         NRECNT = NRECNT + 1
      ELSE IF ( MULTIF.EQ.2 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,0,SUFNRE,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF(NRECNT.EQ.1) THEN
            OPEN(UNIT=CHNREF,STATUS='NEW',FILE=FILBUF)
         ENDIF
         NRECNT = NRECNT + 1
      ENDIF
      GOTO 2
    1 WRITE (CHRSLT,19)
      IF ( MULTIF.EQ.1 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,NRHCNT,SUFNRH,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         OPEN(UNIT=CHNRHF,STATUS='NEW',FILE=FILBUF)
         NRHCNT = NRHCNT + 1
      ELSE IF ( MULTIF.EQ.2 ) THEN
         CALL FILSUF(JOBNAM,FILBUF,0,SUFNRH,IFAIL,DEBUG)
         IF(IFAIL.NE.0) RETURN
         IF(NRHCNT.EQ.1) THEN
            OPEN(UNIT=CHNRHF,STATUS='NEW',FILE=FILBUF)
         ENDIF
         NRHCNT = NRHCNT + 1
      ENDIF

    2 ZNRT = ZNR - DZNR
      DO 11 I = 1 , NRZ
         ZNRT = ZNRT + DZNR
         IF ( NEAR.EQ.0 ) GOTO 3
         CTH = DCOS(TA*ZNRT)
         STH = DSIN(TA*ZNRT)
    3    YNRT = YNR - DYNR
         DO 10 J = 1 , NRY
            YNRT = YNRT + DYNR
            IF ( NEAR.EQ.0 ) GOTO 4
            CPH = DCOS(TA*YNRT)
            SPH = DSIN(TA*YNRT)
    4       XNRT = XNR - DXNR
            DO 9 KK = 1 , NRX
               XNRT = XNRT + DXNR
               IF ( NEAR.EQ.0 ) GOTO 5
               XOB = XNRT*STH*CPH
               YOB = XNRT*STH*SPH
               ZOB = XNRT*CTH
               GOTO 6
    5          XOB = XNRT
               YOB = YNRT
               ZOB = ZNRT
    6          TMP1 = XOB/WLAM
               TMP2 = YOB/WLAM
               TMP3 = ZOB/WLAM
               IF ( NFEH.EQ.1 ) GOTO 7
               CALL NEFLD(TMP1,TMP2,TMP3,EX,EY,EZ,LD,SALP,AIR,AII,
     &                    BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,
     &                    SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,
     &                    DEBUG)
               IF(IFAIL.NE.0) RETURN
               GOTO 8
    7          CALL NHFLD(TMP1,TMP2,TMP3,EX,EY,EZ,LD,SALP,AIR,AII,
     &                    BIR,BII,CIR,CII,CUR,N,M,ICON1,ICON2,X,Y,Z,
     &                    SI,BI,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,IFAIL,
     &                    DEBUG)
               IF(IFAIL.NE.0) RETURN
    8          TMP1 = CDABS(EX)
               TMP2 = CANG(EX)
               TMP3 = CDABS(EY)
               TMP4 = CANG(EY)
               TMP5 = CDABS(EZ)
               TMP6 = CANG(EZ)
 
C     IDF - Peak field calculation from NEC81.
C     WRITE(CHRSLT,11)  XOB,YOB,ZOB,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6
               EZ2 = TMP5*TMP5
               EY2 = TMP3*TMP3
               EX2 = TMP1*TMP1
               AT = EZ2*DCOS(TMP6*2.0D0/RPD) + EY2*DCOS(TMP4*2.0D0/RPD)
     &              + EX2*DCOS(TMP2*2.0D0/RPD)
               BT = EZ2*DSIN(TMP6*2.0D0/RPD) + EY2*DSIN(TMP4*2.0D0/RPD)
     &              + EX2*DSIN(TMP2*2.0D0/RPD)
               CT = AT*AT + BT*BT
               ETOTAL = 0.5D0*(EZ2+EY2+EX2) + 0.5D0*DSQRT(CT)
               ETOTAL = DSQRT(ETOTAL)
C     End of peak field calculation.
C     IDF - Modified to allow output to separate files.
               IF ( MULTIF.GT.0 ) THEN
                  IF( NFEH.EQ.1 ) THEN
                     WRITE (CHNRHF,18) XOB , YOB , ZOB , TMP1 , TMP2 , 
     &                                 TMP3 , TMP4 , TMP5 , TMP6 , 
     &                                 ETOTAL                  
                  ELSE
                     WRITE (CHNREF,18) XOB , YOB , ZOB , TMP1 , TMP2 , 
     &                                 TMP3 , TMP4 , TMP5 , TMP6 , 
     &                                 ETOTAL
                  ENDIF
               ENDIF
C     End IDF.
               WRITE (CHRSLT,18) XOB , YOB , ZOB , TMP1 , TMP2 , 
     &              TMP3 , TMP4 , TMP5 , TMP6 , ETOTAL                  

               IF ( IPLP1.NE.2 ) GOTO 9
               GOTO (14,15,16) , IPLP4
   14          XXX = XOB
               GOTO 17
   15          XXX = YOB
               GOTO 17
   16          XXX = ZOB
   17          CONTINUE
               IF ( IPLP2.NE.2 ) GOTO 13
               IF ( IPLP3.EQ.1 ) WRITE (CHRPAT,*) XXX , TMP1 , TMP2
               IF ( IPLP3.EQ.2 ) WRITE (CHRPAT,*) XXX , TMP3 , TMP4
               IF ( IPLP3.EQ.3 ) WRITE (CHRPAT,*) XXX , TMP5 , TMP6
               IF ( IPLP3.EQ.4 ) WRITE (CHRPAT,*) XXX , TMP1 , TMP2 , 
     &                                    TMP3 , TMP4 , TMP5 , TMP6
C     Added single line for peak field output from NEC81.
               IF ( IPLP3.EQ.5 ) WRITE (CHRPAT,*) XXX , ETOTAL
               GOTO 9
   13          IF ( IPLP2.NE.1 ) GOTO 9
               IF ( IPLP3.EQ.1 ) WRITE (CHRPAT,*) XXX , EX
               IF ( IPLP3.EQ.2 ) WRITE (CHRPAT,*) XXX , EY
               IF ( IPLP3.EQ.3 ) WRITE (CHRPAT,*) XXX , EZ
               IF ( IPLP3.EQ.4 ) WRITE (CHRPAT,*) XXX , EX , EY , EZ
 
    9       CONTINUE
   10    CONTINUE
   11 CONTINUE

      IF ( MULTIF.EQ.1 ) THEN
         IF (NFEH.EQ.1) THEN
            CLOSE(CHNRHF)
         ELSE
            CLOSE(CHNREF)
         ENDIF
      ENDIF

      RETURN
      
C     Format changed for peak field output from NEC81.
   12 FORMAT (///,35X,'- - - NEAR ELECTRIC FIELDS - - -',//,12X,
     &        '-  LOCATION  -',21X,'-  EX  -',15X,'-  EY  -',15X,
     &        '-  EZ  -',/,8X,'X',10X,'Y',10X,'Z',10X,'MAGNITUDE',3X,
     &        'PHASE',6X,'MAGNITUDE',3X,'PHASE',6X,'MAGNITUDE',3X,
     &        'PHASE',/,6X,'METERS',5X,'METERS',5X,'METERS',8X,
     &        'VOLTS/M',3X,'DEGREES',6X,'VOLTS/M',3X,'DEGREES',6X,
     &        'VOLTS/M',3X,'DEGREES',9X,'VOLTS/M')
   18 FORMAT (2X,3(2X,F9.4),1X,3(3X,1P,E11.4,2X,0P,F7.2),6X,E11.4)
   19 FORMAT (///,35X,'- - - NEAR MAGNETIC FIELDS - - -',//,12X,
     &        '-  LOCATION  -',21X,'-  HX  -',15X,'-  HY  -',15X,
     &        '-  HZ  -',/,8X,'X',10X,'Y',10X,'Z',10X,'MAGNITUDE',3X,
     &        'PHASE',6X,'MAGNITUDE',3X,'PHASE',6X,'MAGNITUDE',3X,
     &        'PHASE',/,6X,'METERS',5X,'METERS',5X,'METERS',9X,'AMPS/M',
     &        3X,'DEGREES',7X,'AMPS/M',3X,'DEGREES',7X,'AMPS/M',3X,
     &        'DEGREES',9X,'AMPS/M')
 
      END
