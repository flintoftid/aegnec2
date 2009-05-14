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
 
      SUBROUTINE EFLD( XI , YI , ZI , AI , IJ , IFAIL )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Computes the near E field of a segment with sine, cosine, and      **
C* constant currents in free space or above ground.                   **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XI                                                        **
C* INPUT  - YI                                                        **
C* INPUT  - ZI                                                        **
C* INPUT  - AI                                                        **
C* INPUT  - IJ                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS     ISNOR   SN      XO      XSN     YO     **
C*             YSN     ZO                                             **
C* passes arg  B       IND1    IND2    S                              **
C* uses value  B       CABJ    EXC     EXK     EXS     EYC     EYK    **
C*             EYS     EZC     EZK     EZS     FRATI   IEXK    IND1   **
C*             IND2    IPERF   KSYMP   NRADL   RKH     S       SABJ   **
C*             SALPJ   SCRWL   SN      T1      T2      XJ      XSN    **
C*             YJ      YSN     ZJ      ZRATI                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       EKSC    EKSCX   ROM2    SFLDS                          **
C* called by   CMWW    NEFLD   QDSRC                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL EKSC , EKSCX , ROM2 , SFLDS

C     Dummy arguments.
      INTEGER IJ , IFAIL
      REAL*8 XI , YI , ZI , AI 

C     Local variables.
      INTEGER IIJX , IIP 
      REAL*8 CTH , DMIN , ETA , PI , PX , PY , R , RFL , RH , RHOSPC , 
     &       RHOX , RHOY , RHOZ , RMAG , SALPR , SHAF , TP , XIJ , 
     &       XSPEC , XYMAG , YIJ , YSPEC , ZIJ , ZP
      COMPLEX*16 TXK , TYK , TZK , TXS , TYS , TZS , TXC , TYC , TZC , 
     &           EPX , EPY , REFS , REFPS , ZRSIN , ZRATX , ZSCRN , 
     &           TEZS , TERS , TEZC , TERC , TEZK , TERK
      COMPLEX*16 EGND(9)

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'gnd.inc'
      INCLUDE 'incom.inc'

C     Equivalent local storage.
      EQUIVALENCE (EGND(1),TXK) , (EGND(2),TYK) , (EGND(3),TZK) , 
     &            (EGND(4),TXS) , (EGND(5),TYS) , (EGND(6),TZS) , 
     &            (EGND(7),TXC) , (EGND(8),TYC) , (EGND(9),TZC)
 
C     Data initialisation.
      DATA ETA /376.7303134617706554681984D0/
      DATA PI /3.141592653589793238462643D0/
      DATA TP /6.283185307179586476925286D0/ 
 

      XIJ = XI - XJ
      YIJ = YI - YJ
      IIJX = IJ
      RFL = -1.0D0
      DO 12 IIP = 1 , KSYMP
         IF ( IIP.EQ.2 ) IIJX = 1
         RFL = -RFL
         SALPR = SALPJ*RFL
         ZIJ = ZI - RFL*ZJ
         ZP = XIJ*CABJ + YIJ*SABJ + ZIJ*SALPR
         RHOX = XIJ - CABJ*ZP
         RHOY = YIJ - SABJ*ZP
         RHOZ = ZIJ - SALPR*ZP
         RH = DSQRT(RHOX*RHOX+RHOY*RHOY+RHOZ*RHOZ+AI*AI)
         IF ( RH.GT.1.D-10 ) GOTO 1
         RHOX = 0.0D0
         RHOY = 0.0D0
         RHOZ = 0.0D0
         GOTO 2
    1    RHOX = RHOX/RH
         RHOY = RHOY/RH
         RHOZ = RHOZ/RH
    2    R = DSQRT(ZP*ZP+RH*RH)
         IF ( R.LT.RKH ) GOTO 3

C     Lumped current element approx. for large separations.

         RMAG = TP*R
         CTH = ZP/R
         PX = RH/R
         TXK = DCMPLX(DCOS(RMAG),-DSIN(RMAG))
         PY = TP*R*R
         TYK = ETA*CTH*TXK*DCMPLX(1.0D0,-1.D+0/RMAG)/PY
         TZK = ETA*PX*TXK*DCMPLX(1.0D0,RMAG-1.0D0/RMAG)/(2.0D0*PY)
         TEZK = TYK*CTH - TZK*PX
         TERK = TYK*PX + TZK*CTH
         RMAG = DSIN(PI*S)/PI
         TEZC = TEZK*RMAG
         TERC = TERK*RMAG
         TEZK = TEZK*S
         TERK = TERK*S
         TXS = (0.0D0,0.0D0)
         TYS = (0.0D0,0.0D0)
         TZS = (0.0D0,0.0D0)
         GOTO 6
    3    IF ( IEXK.EQ.1 ) GOTO 4

C     EKSC for thin wire approx. or EKSCX for extended T.W. approx.

         CALL EKSC(S,ZP,RH,TP,IIJX,TEZS,TERS,TEZC,TERC,TEZK,TERK)
         GOTO 5
    4    CALL EKSCX(B,S,ZP,RH,TP,IIJX,IND1,IND2,TEZS,TERS,TEZC,TERC,
     &              TEZK,TERK)
    5    TXS = TEZS*CABJ + TERS*RHOX
         TYS = TEZS*SABJ + TERS*RHOY
         TZS = TEZS*SALPR + TERS*RHOZ
    6    TXK = TEZK*CABJ + TERK*RHOX
         TYK = TEZK*SABJ + TERK*RHOY
         TZK = TEZK*SALPR + TERK*RHOZ
         TXC = TEZC*CABJ + TERC*RHOX
         TYC = TEZC*SABJ + TERC*RHOY
         TZC = TEZC*SALPR + TERC*RHOZ
         IF ( IIP.NE.2 ) GOTO 11
         IF ( IPERF.GT.0 ) GOTO 10
         ZRATX = ZRATI
         RMAG = R
         XYMAG = DSQRT(XIJ*XIJ+YIJ*YIJ)

C     Set parameters for radial wire ground screen.

         IF ( NRADL.EQ.0 ) GOTO 7
         XSPEC = (XI*ZJ+ZI*XJ)/(ZI+ZJ)
         YSPEC = (YI*ZJ+ZI*YJ)/(ZI+ZJ)
         RHOSPC = DSQRT(XSPEC*XSPEC+YSPEC*YSPEC+T2*T2)
         IF ( RHOSPC.GT.SCRWL ) GOTO 7
         ZSCRN = T1*RHOSPC*DLOG(RHOSPC/T2)
         ZRATX = (ZSCRN*ZRATI)/(ETA*ZRATI+ZSCRN)
    7    IF ( XYMAG.GT.1.D-6 ) GOTO 8

C     Calculation of reflection coefficients when ground is specified.

         PX = 0.0D0
         PY = 0.0D0
         CTH = 1.0D0
         ZRSIN = (1.0D0,0.0D0)
         GOTO 9
    8    PX = -YIJ/XYMAG
         PY = XIJ/XYMAG
         CTH = ZIJ/RMAG
         ZRSIN = CDSQRT(1.0D0-ZRATX*ZRATX*(1.0D0-CTH*CTH))
    9    REFS = (CTH-ZRATX*ZRSIN)/(CTH+ZRATX*ZRSIN)
         REFPS = -(ZRATX*CTH-ZRSIN)/(ZRATX*CTH+ZRSIN)
         REFPS = REFPS - REFS
         EPY = PX*TXK + PY*TYK
         EPX = PX*EPY
         EPY = PY*EPY
         TXK = REFS*TXK + REFPS*EPX
         TYK = REFS*TYK + REFPS*EPY
         TZK = REFS*TZK
         EPY = PX*TXS + PY*TYS
         EPX = PX*EPY
         EPY = PY*EPY
         TXS = REFS*TXS + REFPS*EPX
         TYS = REFS*TYS + REFPS*EPY
         TZS = REFS*TZS
         EPY = PX*TXC + PY*TYC
         EPX = PX*EPY
         EPY = PY*EPY
         TXC = REFS*TXC + REFPS*EPX
         TYC = REFS*TYC + REFPS*EPY
         TZC = REFS*TZC
   10    EXK = EXK - TXK*FRATI
         EYK = EYK - TYK*FRATI
         EZK = EZK - TZK*FRATI
         EXS = EXS - TXS*FRATI
         EYS = EYS - TYS*FRATI
         EZS = EZS - TZS*FRATI
         EXC = EXC - TXC*FRATI
         EYC = EYC - TYC*FRATI
         EZC = EZC - TZC*FRATI
         GOTO 12
   11    EXK = TXK
         EYK = TYK
         EZK = TZK
         EXS = TXS
         EYS = TYS
         EZS = TZS
         EXC = TXC
         EYC = TYC
         EZC = TZC
   12 CONTINUE
      IF ( IPERF.EQ.2 ) GOTO 13
      RETURN

C     Field due to ground using Sommerfeld/Norton.

   13 SN = DSQRT(CABJ*CABJ+SABJ*SABJ)
      IF ( SN.LT.1.D-5 ) GOTO 14
      XSN = CABJ/SN
      YSN = SABJ/SN
      GOTO 15
   14 SN = 0.0D0
      XSN = 1.0D0
      YSN = 0.0D0

C     Displace observation point for thin wire approximation.

   15 ZIJ = ZI + ZJ
      SALPR = -SALPJ
      RHOX = SABJ*ZIJ - SALPR*YIJ
      RHOY = SALPR*XIJ - CABJ*ZIJ
      RHOZ = CABJ*YIJ - SABJ*XIJ
      RH = RHOX*RHOX + RHOY*RHOY + RHOZ*RHOZ
      IF ( RH.GT.1.D-10 ) GOTO 16
      XO = XI - AI*YSN
      YO = YI + AI*XSN
      ZO = ZI
      GOTO 17
   16 RH = AI/DSQRT(RH)
      IF ( RHOZ.LT.0. ) RH = -RH
      XO = XI + RH*RHOX
      YO = YI + RH*RHOY
      ZO = ZI + RH*RHOZ
   17 R = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ
      IF ( R.GT.0.95D0 ) GOTO 18

C     Field from interpolation is integrated over segment.

      ISNOR = 1
      DMIN = DBLE(EXK*DCONJG(EXK) + EYK*DCONJG(EYK) + EZK*DCONJG(EZK))
      DMIN = 0.01D0*DSQRT(DMIN)
      SHAF = 0.5D0*S
      CALL ROM2(-SHAF,SHAF,EGND,DMIN,IFAIL)
      IF(IFAIL.NE.0) RETURN
      GOTO 19

C     Norton field equations and lumped current element approximation.

   18 ISNOR = 2
      CALL SFLDS(0.0D0,EGND)
      GOTO 22
   19 ZP = XIJ*CABJ + YIJ*SABJ + ZIJ*SALPR
      RH = R - ZP*ZP
      IF ( RH.GT.1.D-10 ) GOTO 20
      DMIN = 0.0D0
      GOTO 21
   20 DMIN = DSQRT(RH/(RH+AI*AI))
   21 IF ( DMIN.GT.0.95D0 ) GOTO 22
      PX = 1.0D0 - DMIN
      TERK = (TXK*CABJ+TYK*SABJ+TZK*SALPR)*PX
      TXK = DMIN*TXK + TERK*CABJ
      TYK = DMIN*TYK + TERK*SABJ
      TZK = DMIN*TZK + TERK*SALPR
      TERS = (TXS*CABJ+TYS*SABJ+TZS*SALPR)*PX
      TXS = DMIN*TXS + TERS*CABJ
      TYS = DMIN*TYS + TERS*SABJ
      TZS = DMIN*TZS + TERS*SALPR
      TERC = (TXC*CABJ+TYC*SABJ+TZC*SALPR)*PX
      TXC = DMIN*TXC + TERC*CABJ
      TYC = DMIN*TYC + TERC*SABJ
      TZC = DMIN*TZC + TERC*SALPR
   22 EXK = EXK + TXK
      EYK = EYK + TYK
      EZK = EZK + TZK
      EXS = EXS + TXS
      EYS = EYS + TYS
      EZS = EZS + TZS
      EXC = EXC + TXC
      EYC = EYC + TYC
      EZC = EZC + TZC
 
      RETURN
 
      END
