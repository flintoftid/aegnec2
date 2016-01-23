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
C* ---------------------------------------------------------------------
C*
C* Contains I/O channel numbers, file name suffices and other parameters
C* for use in the FORTRAN code.
C*

C     Patch offset in ICON arrays.
      INTEGER PATOFF
      PARAMETER( PATOFF = 10000 )

C     Input file suffix. This must match the definition of SUFINP 
C     in nec.h.
      CHARACTER*3 SUFINP
      PARAMETER( SUFINP = 'nec' )

C     Output file suffix.
      CHARACTER*3 SUFRES
      PARAMETER( SUFRES = 'res' )

C     Radiation pattern file suffix.
      CHARACTER*3 SUFRDP
      PARAMETER( SUFRDP = 'rdp' )

C     Near F field file suffix.
      CHARACTER*3 SUFNRE
      PARAMETER( SUFNRE = 'nre' )

C     Near H field file suffix.
      CHARACTER*3 SUFNRH
      PARAMETER( SUFNRH = 'nrh' )

C     Impedance file suffix.
      CHARACTER*3 SUFIMP
      PARAMETER( SUFIMP = 'imp' )

C     NGF file suffix.
      CHARACTER*3 SUFNGF
      PARAMETER( SUFNGF = 'ngf' )

C     Sommerfeld ground file suffix.
      CHARACTER*3 SUFSOM
      PARAMETER( SUFSOM = 'som' )

C     Input file channel.
      INTEGER CHINPT
      PARAMETER( CHINPT = 2 )

C     Results file channel.
      INTEGER CHRSLT
      PARAMETER( CHRSLT = 3 )

C     Radiation pattern files.
      INTEGER CHRPAT
      PARAMETER( CHRPAT = 8 )

C     Temporary file channels.
      INTEGER CHTMP1 , CHTMP2 , CHTMP3 , CHTMP4 , CHTMP5 , CHTMP6
      PARAMETER( CHTMP1 = 11 , CHTMP2 = 12 , CHTMP3 = 13 ,
     &           CHTMP4 = 14 , CHTMP5 = 15 , CHTMP6 = 16 )

C     NGF file channel.
      INTEGER CHNGFL
      PARAMETER( CHNGFL = 20 )

C     SOMNEC file channel.
      INTEGER CHSOMN
      PARAMETER( CHSOMN = 21 )

C     Near E field file.
      INTEGER CHNREF
      PARAMETER( CHNREF = 81 )

C     Near H field file.
      INTEGER CHNRHF
      PARAMETER( CHNRHF = 82 )

C     Radiation pattern file.
      INTEGER CHRDPT
      PARAMETER( CHRDPT = 45 )

C     Impedance file.
      INTEGER CHIMPD
      PARAMETER( CHIMPD = 78 )

