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
C* EVLCOM common block.                                               **
C*                                                                    **
C*--------------------------------------------------------------------**

      INTEGER JH
      REAL*8 CK2 , CK2SQ , TKMAG , TSMAG , CK1R , ZPH , RHO  
      COMPLEX*16 CKSM , CT1 , CT2 , CT3 , CK1 , CK1SQ

      COMMON /EVLCOM/ CKSM , CT1 , CT2 , CT3 , CK1 , CK1SQ , CK2 , 
     &                CK2SQ , TKMAG , TSMAG , CK1R , ZPH , RHO , JH

