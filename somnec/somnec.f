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

      BLOCK DATA

C*--------------------------------------------------------------------**
C*                                                                    **
C* Initialise GGRID common block data.                                **
C*                                                                    **
C*--------------------------------------------------------------------**

      INCLUDE 'ggrid.inc'

      DATA NXA /11 , 17 , 9/
      DATA NYA /10 , 5 , 8/
      DATA XSA /0.0D0 , 0.2D0 , 0.2D0/
      DATA YSA /0.0D0 , 0.0D0 , 0.3490658504D0/
      DATA DXA /0.02D0 , 0.05D0 , 0.1D0/
      DATA DYA /0.1745329252D0 , 0.0872654626D0 , 0.1745329252D0/

      END

C*--------------------------------------------------------------------**

      SUBROUTINE SOMNEC( IFAIL , VERBSE , DEBUG , IPT , EPR , SIG , 
     &                   FMHZ , OTFILE )

C*--------------------------------------------------------------------**
C*                                                                    **
C* PROGRAM TO GENERATE NEC INTERPOLATION GRIDS FOR FIELDS DUE TO      **
C* GROUND.  FIELD COMPONENTS ARE COMPUTED BY NUMERICAL EVALUATION     **
C* OF MODIFIED SOMMERFELD INTEGRALS.                                  **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - IFAIL                                                     **
C* INPUT - VERBSE                                                     **
C* INPUT - DEBUG                                                      **
C* INPUT - EPR                                                        **
C* INPUT - SIG                                                        **
C* INPUT - FMHZ                                                       **
C* INPUT - OTFILE                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    AR1     AR2     AR3     CK1     CK1R    CK1SQ   CK2    **
C*             CK2SQ   CKSM    CT1     CT2     CT3     DXA     DYA    **
C*             EPSCF   NXA     NYA     RHO     TKMAG   TSMAG   XSA    **
C*             YSA     ZPH                                            **
C* uses value  /GGRID/ AR1     AR2     AR3     CK1     CK1SQ   CK2    **
C*             CK2SQ   DXA     DYA     EPSCF   NXA     NYA     RHO    **
C*             XSA     YSA     ZPH                                    **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       EVLUA   TIMER                                          **
C* called by   ** NOTHING **                                          **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL EVLUA , TIMER

C     Dummy arguments.
      CHARACTER*32 OTFILE
      INTEGER IFAIL , VERBSE , DEBUG , IPT
      REAL*8 EPR , FMHZ , SIG 

C     Local variables.
      CHARACTER*3  LCOMP(4)
      INTEGER IR , IRS , ITH , K , L , NR , NTH 
      REAL*8 DR , DTH , R , RK , TFAC1 , TFAC2 , THET , WLAM , TIM, 
     &       TST, TFN 
      COMPLEX*16  ERV , EZV , ERH , EPH , CL1 , CL2 , CON

C     Common storage.
      INCLUDE 'evlcom.inc'
      INCLUDE 'ggrid.inc'

C     Data initialisation.
      DATA LCOMP /'ERV' , 'EZV' , 'ERH' , 'EPH'/
 
      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING SOMNEC'
      ENDIF

      IF (SIG.LT.0) GO TO 1
      WLAM=299.8D0/FMHZ
      EPSCF=DCMPLX(EPR,-SIG*WLAM*59.96D0)
      GO TO 2
1     EPSCF=DCMPLX(EPR,SIG)
C2    TST=SECNDS(0.0)
2     CALL TIMER(TST)
      CK2=6.283185308D0
      CK2SQ=CK2*CK2

C     Sommerfeld integral evaluation uses EXP(-JWT), NEC uses EXP(+JWT),
C     hence need DCONJG(EPSCF).  Conjugate of fields occurs in subroutine
C     EVALUA.

      CK1SQ=CK2SQ*DCONJG(EPSCF)
      CK1=CDSQRT(CK1SQ)
      CK1R=DBLE(CK1)
      TKMAG=100.0D0*CDABS(CK1)
      TSMAG=DBLE(100.0D0*CK1*DCONJG(CK1))
      CKSM=CK2SQ/(CK1SQ+CK2SQ)
      CT1=0.5D0*(CK1SQ-CK2SQ)
      ERV=CK1SQ*CK1SQ
      EZV=CK2SQ*CK2SQ
      CT2=0.125D0*(ERV-EZV)
      ERV=ERV*CK1SQ
      EZV=EZV*CK2SQ
      CT3=.0625D0*(ERV-EZV)

C     Loop over 3 grid regions.

      DO 6 K=1,3
      NR=NXA(K)
      NTH=NYA(K)
      DR=DXA(K)
      DTH=DYA(K)
      R=XSA(K)-DR
      IRS=1
      IF (K.EQ.1) R=XSA(K)
      IF (K.EQ.1) IRS=2

C     Loop over R (R=SQRT(RHO**2 + (Z+H)**2)).

      DO 6 IR=IRS,NR
      R=R+DR
      THET=YSA(K)-DTH

C     Loop over theta (THETA=ATAN((Z+H)/RHO)).

      DO 6 ITH=1,NTH
      THET=THET+DTH
      RHO=R*DCOS(THET)
      ZPH=R*DSIN(THET)
      IF (RHO.LT.1.0D-7) RHO=1.0D-8
      IF (ZPH.LT.1.0D-7) ZPH=0.0D0
      CALL EVLUA (ERV,EZV,ERH,EPH)
      RK=CK2*R
      CON=-(0.0D0,4.77147D0)*R/DCMPLX(DCOS(RK),-DSIN(RK))
      GO TO (3,4,5), K
3     AR1(IR,ITH,1)=ERV*CON
      AR1(IR,ITH,2)=EZV*CON
      AR1(IR,ITH,3)=ERH*CON
      AR1(IR,ITH,4)=EPH*CON
      GO TO 6
4     AR2(IR,ITH,1)=ERV*CON
      AR2(IR,ITH,2)=EZV*CON
      AR2(IR,ITH,3)=ERH*CON
      AR2(IR,ITH,4)=EPH*CON
      GO TO 6
5     AR3(IR,ITH,1)=ERV*CON
      AR3(IR,ITH,2)=EZV*CON
      AR3(IR,ITH,3)=ERH*CON
      AR3(IR,ITH,4)=EPH*CON
6     CONTINUE

C     Fill grid 1 for R equal to zero.

      CL2=-(0.0D0,188.37D0)*(EPSCF-1.0D0)/(EPSCF+1.0D0)
      CL1=CL2/(EPSCF+1.0D0)
      EZV=EPSCF*CL1
      THET=-DTH
      NTH=NYA(1)
      DO 9 ITH=1,NTH
      THET=THET+DTH
      IF (ITH.EQ.NTH) GO TO 7
      TFAC2=DCOS(THET)
      TFAC1=(1.0D0-DSIN(THET))/TFAC2
      TFAC2=TFAC1/TFAC2
      ERV=EPSCF*CL1*TFAC1
      ERH=CL1*(TFAC2-1.0D0)+CL2
      EPH=CL1*TFAC2-CL2
      GO TO 8
7     ERV=0.0D0
      ERH=CL2-0.5D0*CL1
      EPH=-ERH
8     AR1(1,ITH,1)=ERV
      AR1(1,ITH,2)=EZV
      AR1(1,ITH,3)=ERH
9     AR1(1,ITH,4)=EPH
C     TIM=SECNDS (TST)
      CALL TIMER(TFN)
      TIM = TFN-TST

C     Write grid on TAPE21.
      IF( OTFILE(1:5).NE.'NOOUT' ) THEN
         OPEN (UNIT=21,FILE=OTFILE,FORM='UNFORMATTED',STATUS='NEW',
     &         ERR=21)
         WRITE (21) AR1,AR2,AR3,EPSCF,DXA,DYA,XSA,YSA,NXA,NYA
         CLOSE (UNIT=21)
      ENDIF
      IF (IPT.EQ.0) GO TO 14

C     Print grid.

      PRINT 17, EPSCF
      DO 13 K=1,3
      NR=NXA(K)
      NTH=NYA(K)
      PRINT 18, K,XSA(K),DXA(K),NR,YSA(K),DYA(K),NTH
      DO 13 L=1,4
      PRINT 19, LCOMP(L)
      DO 13 IR=1,NR
      GO TO (10,11,12), K
10    PRINT 20, IR,(AR1(IR,ITH,L),ITH=1,NTH)
      GO TO 13
11    PRINT 20, IR,(AR2(IR,ITH,L),ITH=1,NTH)
      GO TO 13
12    PRINT 20, IR,(AR3(IR,ITH,L),ITH=1,NTH)
13    CONTINUE
14    CONTINUE
      IF( VERBSE.EQ.1 ) THEN
      PRINT 16, TIM
      ENDIF
      GO TO 23
21    IFAIL=1
23    RETURN

 16   FORMAT (' TIME=',E12.3,' SECONDS')
 17   FORMAT (' NEC GROUND INTERPOLATION GRID',/,
     &        ' DIELECTRIC CONSTANT=',2E12.5)
 18   FORMAT (///,' GRID',I2,/,4X,'R(1)=',F7.4,4X,'DR=',F7.4,4X,'NR=',
     &        I3,/,' THET(1)=',F7.4,3X,'DTH=',F7.4,3X,'NTH=',I3,//)
 19   FORMAT (///1X,A3)
 20   FORMAT (' IR=',I3,/1X,(10(1PE12.5)))
 24   FORMAT (A)

      END
