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
C* /GGRID/ common block.                                              **   
C*                                                                    **
C*----------------------Aggregated COMMON block usage-----------------**
C*                                                                    **
C* /GGRID/ - changed by GFIL     NEC2D                                **
C*           EQV/arg in INTRP                                         **
C*              used by GFOUT    INTRP    NEC2D                       **
C*                                                                    **
C*--------------------------------------------------------------------**

      INTEGER NXA , NYA
      REAL*8 DXA , DYA , XSA , YSA
      COMPLEX*16 AR1 , AR2 , AR3 , EPSCF

      COMMON /GGRID / AR1(11,10,4) , AR2(17,5,4) , AR3(9,8,4) , EPSCF , 
     &                DXA(3) , DYA(3) , XSA(3) , YSA(3) , NXA(3) , 
     &                NYA(3)
