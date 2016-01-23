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

      SUBROUTINE HSFLD( XI , YI , ZI , AI )

C*--------------------------------------------------------------------**
C*                                                                    **
C* HSFLD computes the H field for constant, sine, and cosine current  **
C* on a segment including ground effects.                             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XI                                                        **
C* INPUT  - YI                                                        **
C* INPUT  - ZI                                                        **
C* INPUT  - AI                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    EXC     EXK     EXS     EYC     EYK     EYS     EZC    **
C*             EZK     EZS                                            **
C* passes arg  S                                                      **
C* uses value  CABJ    EXC     EXK     EXS     EYC     EYK     EYS    **
C*             EZC     EZK     EZS     IPERF   KSYMP   NRADL   S      **
C*             SABJ    SALPJ   SCRWL   T1      T2      XJ      YJ     **
C*             ZJ      ZRATI                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       HSFLX                                                  **
C* called by   CMWS    NHFLD   QDSRC                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL HSFLX

C     Dummy arguments.
      REAL*8 XI , YI , ZI , AI

C     Local variables.
      INTEGER IIP
      REAL*8 CTH , ETA , PHX , PHY , PHZ , PX , PY , RFL , RH , 
     &       RHOSPC , RHOX , RHOY , RHOZ , RMAG , SALPR , XIJ , 
     &       XSPEC , XYMAG , YIJ , YSPEC , ZIJ , ZP
      COMPLEX*16 HPK , HPS , HPC , QX , QY , QZ , RRV , RRH , ZRATX

C     Common storage.
      INCLUDE 'dataj.inc'
      INCLUDE 'gnd.inc'
 
C     Data initialisation.
      DATA ETA /376.7303134617706554681984D0/
 

      XIJ = XI - XJ
      YIJ = YI - YJ
      RFL = -1.0D0
      DO 7 IIP = 1 , KSYMP
         RFL = -RFL
         SALPR = SALPJ*RFL
         ZIJ = ZI - RFL*ZJ
         ZP = XIJ*CABJ + YIJ*SABJ + ZIJ*SALPR
         RHOX = XIJ - CABJ*ZP
         RHOY = YIJ - SABJ*ZP
         RHOZ = ZIJ - SALPR*ZP
         RH = DSQRT(RHOX*RHOX+RHOY*RHOY+RHOZ*RHOZ+AI*AI)
         IF ( RH.GT.1.D-10 ) GOTO 1
         EXK = 0.0D0
         EYK = 0.0D0
         EZK = 0.0D0
         EXS = 0.0D0
         EYS = 0.0D0
         EZS = 0.0D0
         EXC = 0.0D0
         EYC = 0.0D0
         EZC = 0.0D0
         GOTO 7
    1    RHOX = RHOX/RH
         RHOY = RHOY/RH
         RHOZ = RHOZ/RH
         PHX = SABJ*RHOZ - SALPR*RHOY
         PHY = SALPR*RHOX - CABJ*RHOZ
         PHZ = CABJ*RHOY - SABJ*RHOX
         CALL HSFLX(S,RH,ZP,HPK,HPS,HPC)
         IF ( IIP.NE.2 ) GOTO 6
         IF ( IPERF.EQ.1 ) GOTO 5
         ZRATX = ZRATI
         RMAG = DSQRT(ZP*ZP+RH*RH)
         XYMAG = DSQRT(XIJ*XIJ+YIJ*YIJ)

C     Set parameters for radial wire ground screen.

         IF ( NRADL.EQ.0 ) GOTO 2
         XSPEC = (XI*ZJ+ZI*XJ)/(ZI+ZJ)
         YSPEC = (YI*ZJ+ZI*YJ)/(ZI+ZJ)
         RHOSPC = DSQRT(XSPEC*XSPEC+YSPEC*YSPEC+T2*T2)
         IF ( RHOSPC.GT.SCRWL ) GOTO 2
         RRV = T1*RHOSPC*DLOG(RHOSPC/T2)
         ZRATX = (RRV*ZRATI)/(ETA*ZRATI+RRV)
    2    IF ( XYMAG.GT.1.D-6 ) GOTO 3

C     Calculation of reflection coefficients when ground is specified.

         PX = 0.0D0
         PY = 0.0D0
         CTH = 1.0D0
         RRV = (1.0D0,0.0D0)
         GOTO 4
    3    PX = -YIJ/XYMAG
         PY = XIJ/XYMAG
         CTH = ZIJ/RMAG
         RRV = CDSQRT(1.0D0-ZRATX*ZRATX*(1.0D0-CTH*CTH))
    4    RRH = ZRATX*CTH
         RRH = -(RRH-RRV)/(RRH+RRV)
         RRV = ZRATX*RRV
         RRV = (CTH-RRV)/(CTH+RRV)
         QY = (PHX*PX+PHY*PY)*(RRV-RRH)
         QX = QY*PX + PHX*RRH
         QY = QY*PY + PHY*RRH
         QZ = PHZ*RRH
         EXK = EXK - HPK*QX
         EYK = EYK - HPK*QY
         EZK = EZK - HPK*QZ
         EXS = EXS - HPS*QX
         EYS = EYS - HPS*QY
         EZS = EZS - HPS*QZ
         EXC = EXC - HPC*QX
         EYC = EYC - HPC*QY
         EZC = EZC - HPC*QZ
         GOTO 7
    5    EXK = EXK - HPK*PHX
         EYK = EYK - HPK*PHY
         EZK = EZK - HPK*PHZ
         EXS = EXS - HPS*PHX
         EYS = EYS - HPS*PHY
         EZS = EZS - HPS*PHZ
         EXC = EXC - HPC*PHX
         EYC = EYC - HPC*PHY
         EZC = EZC - HPC*PHZ
         GOTO 7
    6    EXK = HPK*PHX
         EYK = HPK*PHY
         EZK = HPK*PHZ
         EXS = HPS*PHX
         EYS = HPS*PHY
         EZS = HPS*PHZ
         EXC = HPC*PHX
         EYC = HPC*PHY
         EZC = HPC*PHZ
    7 CONTINUE
 
      RETURN
 
      END
