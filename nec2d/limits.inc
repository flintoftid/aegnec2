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
C* limits.inc - Array limits for static versions of NEC2.             **
C*                                                                    **
C*--------------------------------------------------------------------**

C     Maximum number of segments or patches. Program will accept up
C     to LD segments OR LD patches. Mixture of patches and segments 
C     limited by N+M = LD, N+2M = 2*LD and N+3*M = 3* LD, where 
C     N = number of segments and M = number of patches.
      INTEGER LD
      PARAMETER( LD = 2000 )

C     In core matrix storage. IRESRV is the number of words of core
C     storage for the interaction matrix. Matrix will fit in core if
C     (N+2*M)*(NP+2*MP) <= IRESRV. For out of core solution IRESRV
C     must be at least 2*(N+2*M). For guaranteed in core storage
C     require IRESRV = 4*LD^2.
      INTEGER IRESRV
C     Use this for gguaranteed in core solution with patches & segments:
C     PARAMETER ( IRESRV = 4 * LD * LD )
C     Use this for guaranteed in core solution with only segments:
      PARAMETER ( IRESRV = LD * LD )
C     Or set it manually:
C     PARAMETER ( IRESRV = 4000000 )

C     Maximum number of voltage sources. A model may use up to
C     NSMAX applied field sources and NSMAX current-slope 
C     discontinuity sources.
      INTEGER NSMAX
      PARAMETER( NSMAX = 30 )

C     Maximum number of loading cards.
      INTEGER LOADMX
      PARAMETER( LOADMX = 30 )

C     Maximium number of non-radiating networks.
      INTEGER NETMX
      PARAMETER( NETMX = 30 )

C     Maximum number of frequencies for normalised impedance or
C     maximum number of angles for which received signal strength is
C     stored.
      INTEGER NORMF
      PARAMETER( NORMF = 200 )

C     Maximum number of points in coupling calcuation.
      INTEGER MXCOUP , MXCOU2
      PARAMETER( MXCOUP = 5 )
      PARAMETER( MXCOU2 = MXCOUP * MXCOUP - MXCOUP )

C     Maximum number of segments joining at junctions.
      INTEGER JMAX
      PARAMETER( JMAX = 30 )

C     Maximum number of NGF segments to which new segments or
C     patches connect.
      INTEGER NSMAXX
      PARAMETER( NSMAXX = 50 )

C     Maximum number of NGF patches to which new segments connect.
      INTEGER NPMAX
      PARAMETER( NPMAX = 10 )

