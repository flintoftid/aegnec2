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
C*                                                                    **
C* /MATPAR/ common block.                                             **   
C*                                                                    **
C*----------------------Aggregated COMMON block usage-----------------**
C*                                                                    **
C* /MATPAR/ - changed by FACGF    FBLOCK   FBNGF    GFIL     NEC2D    **
C*                       SOLGF                                        **
C*               used by CMNGF    CMSET    FACGF    FACIO    FACTRS   **
C*                       FBLOCK   FBNGF    GFIL     GFOUT    LFACTR   **
C*                       LTSOLV   LUNSCR   NEC2D    REBLK    SOLGF    **
C*                       SOLVES                                       **
C*                                                                    **
C*--------------------------------------------------------------------**

      INTEGER ICASE , NBLOKS , NPBLK , NLAST , NBLSYM , NPSYM , NLSYM , 
     &        IMAT , ICASX , NBBX , NPBX , NLBX , NBBL , NPBL , NLBL

      COMMON /MATPAR/ ICASE , NBLOKS , NPBLK , NLAST , NBLSYM , NPSYM , 
     &                NLSYM , IMAT , ICASX , NBBX , NPBX , NLBX , NBBL , 
     &                NPBL , NLBL

