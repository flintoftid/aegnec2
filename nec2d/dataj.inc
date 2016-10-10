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
C*                                                                    **
C* /DATAJ/ common block.                                              **   
C*                                                                    **
C*----------------------Aggregated COMMON block usage-----------------**
C*                                                                    **
C* /DATAJ/ - changed by CMNGF    CMSET    CMSS     CMSW     CMWS      **
C*                      CMWW     EFLD     HINTG    HSFLD    NEFLD     **
C*                      NHFLD    PCINT    QDSRC    UNERE              **
C*              used by CMNGF    CMSS     CMSW     CMWS     CMWW      **
C*                      EFLD     HINTG    HSFLD    NEFLD    NHFLD     **
C*                      PCINT    QDSRC    SFLDS    UNERE              **
C*                                                                    **
C*--------------------------------------------------------------------**

      INTEGER IEXK , IND1 , IND2 , IPGND 
      REAL*8 S , B , XJ , YJ , ZJ , CABJ , SABJ , SALPJ , RKH
      COMPLEX *16 EXK , EYK , EZK , EXS , EYS , EZS , EXC , EYC , EZC 

      COMMON /DATAJ/ EXK , EYK , EZK , EXS , EYS , EZS , EXC , EYC , 
     &               EZC , S , B , XJ , YJ , ZJ , CABJ , SABJ , 
     &               SALPJ , RKH , IND1 , IND2 , IPGND , IEXK  

