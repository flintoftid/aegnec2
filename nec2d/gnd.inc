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
C*                                                                    **
C* /GND/ common block.                                                **   
C*                                                                    **
C*----------------------Aggregated COMMON block usage-----------------**
C*                                                                    **
C* /GND/ - changed by GFIL     NEC2D    RDPAT                         **
C*            used by CMSW     EFLD     ETMNS    FFLD     GFIL        **
C*                             GFOUT    HINTG    HSFLD    NEC2D       **
C*                             NEFLD    NHFLD    RDPAT    SFLDS       **
C*                             UNERE                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      INTEGER NRADL , KSYMP , IFAR , IPERF
      REAL*8 CL , CH , SCRWL , SCRWR , T2
      COMPLEX*16 ZRATI , ZRATI2 , FRATI , T1

      COMMON /GND/ ZRATI , ZRATI2 , FRATI , T1 , CL , CH , SCRWL , 
     &             SCRWR , T2 , NRADL , KSYMP , IFAR , IPERF
