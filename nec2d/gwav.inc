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
C* /GWAV/ common block.                                               **   
C*                                                                    **
C*----------------------Aggregated COMMON block usage-----------------**
C*                                                                    **
C* /GWAV/ - changed by GFLD     NEC2D    SFLDS                        **
C*             used by GFLD     GWAVE    NEC2D    SFLDS               **
C*                                                                    **
C*--------------------------------------------------------------------**

      REAL*8 R1 , R2 , ZMH , ZPH
      COMPLEX*16 U , U2 , XX1 , XX2 

      COMMON /GWAV/ U , U2 , XX1 , XX2 , R1 , R2 , ZMH , ZPH

