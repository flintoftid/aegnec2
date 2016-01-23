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

      INTEGER FUNCTION ISEGNO( ITAGI , MX , LD , ITAG, N , IFAIL , 
     &                         DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* ISEGNO returns the segment number of the Mth segment having the    **
C* tag number ITAGI. If ITAGI=0 segment number M is returned.         **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - ITAGI                                                     **
C* INPUT  - MX                                                        **
C* INPUT  - LD                                                        **
C* INPUT  - ITAG                                                      **
C* INPUT  - N                                                         **
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
C* called by   COUPLE  MOVE    NEC2D                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER ITAGI , MX , LD , N , IFAIL , DEBUG
      INTEGER ITAG(2*LD)

C     Local variables.
      INTEGER I , ICNT 


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING ISEGNO'
      ENDIF
 
      IF ( MX.GT.0 ) GOTO 1
      WRITE (CHRSLT,6)
      ISEGNO=0
      IFAIL=22
      RETURN
    1 ICNT = 0
      IF ( ITAGI.NE.0 ) GOTO 2
      ISEGNO = MX
      RETURN
    2 IF ( N.LT.1 ) GOTO 4
      DO 3 I = 1 , N
         IF ( ITAG(I).NE.ITAGI ) GOTO 3
         ICNT = ICNT + 1
         IF ( ICNT.EQ.MX ) GOTO 5
    3 CONTINUE
    4 WRITE (CHRSLT,7) ITAGI
      ISEGNO=0
      IFAIL=23
      RETURN
    5 ISEGNO = I
 
      RETURN
 
    6 FORMAT (4X,
     &'CHECK DATA, PARAMETER SPECIFYING SEGMENT POSITION IN A GROUP OF E
     &QUAL TAGS MUST NOT BE ZERO')
    7 FORMAT (///,10X,'NO SEGMENT HAS AN ITAG OF ',I5)
 
      END
