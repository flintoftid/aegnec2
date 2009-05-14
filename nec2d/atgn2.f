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
 
      REAL*8 FUNCTION ATGN2 ( XX , YY )
 
C*--------------------------------------------------------------------**
C*                                                                    **
C* ATGN2 is arctangent function modified to return 0 when XX=YY=0.    **
C* Most standard ATAN2 functions give an error when XX=YY=0.          **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - XX                                                        **
C* INPUT  - YY                                                        **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   CANG    DATAGN  RDPAT                                  **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     Dummy arguments.
      REAL*8 XX , YY


      IF ( XX ) 3 , 1 , 3
    1 IF ( YY ) 3 , 2 , 3
    2 ATGN2 = 0.0D0
       RETURN
    3 ATGN2 = DATAN2(XX,YY)
 
      RETURN

      END
