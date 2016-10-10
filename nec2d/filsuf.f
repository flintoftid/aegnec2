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

      SUBROUTINE FILSUF( JOBNAM , FILBUF , NUMBER , SUFFIX , IFAIL , 
     &                   DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Construct file name from jobname, file number and file suffix.     **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - JOBNAM                                                    **
C* OUTPUT - FILBUF                                                    **
C* INPUT  - NUMBER                                                    **
C* INPUT  - SUFFIX                                                    **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FDRIVE  GFIL    GFOUT   NEC2D   NFPAT   RDPAT          **
C*                                                                    **
C*--------------------------------------------------------------------**
      
      IMPLICIT NONE

C     Dummy arguments.
      CHARACTER*(*) JOBNAM , FILBUF , SUFFIX
      INTEGER NUMBER , IFAIL , DEBUG 
      
C     Local variables.
      CHARACTER*1 DIGIT(0:9)
      INTEGER D1 , D2 , D3 , D4 , ISPACE

C     Data initialisation
      DATA DIGIT /'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING FILSUF'
      ENDIF

      ISPACE = INDEX( JOBNAM, ' ' )
      IF( ISPACE.EQ.0 ) THEN
         ISPACE=LEN(JOBNAM)
      ELSE
         ISPACE=ISPACE-2
      ENDIF
      
      IF( (LEN(JOBNAM)+LEN(SUFFIX)+4+2).GT.LEN(FILBUF) ) THEN
         IFAIL=19
         RETURN
      ELSE IF( NUMBER.GT.9999 ) THEN
         IFAIL=20
         RETURN
      ELSE IF( NUMBER.LT.0 ) THEN
         IFAIL=21
         RETURN
      ELSE IF( NUMBER.EQ.0 ) THEN
         FILBUF = JOBNAM(1:ISPACE) // '.' // SUFFIX
         RETURN
      ELSE
         D1 = NUMBER / 1000
         D2 = ( NUMBER - D1 * 1000 ) / 100 
         D3 = ( NUMBER - D1 * 1000 - D2 * 100 ) / 10 
         D4 = ( NUMBER - D1 * 1000 - D2 * 100 - D3 * 10 )
         FILBUF = JOBNAM(1:ISPACE) // '-' // DIGIT(D1) // DIGIT(D2) 
     &                             // DIGIT(D3) // DIGIT(D4) // '.' 
     &                             // SUFFIX
      ENDIF

      RETURN

      END

