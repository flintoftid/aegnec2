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

      SUBROUTINE READMN( INUNIT , CODE , I1 , I2 , I3 , I4 , F1 , F2 ,
     &                   F3 , F4 , F5 , F6 , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* READMN reads a control record and parses it.                       **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - INUNIT                                                    **
C* OUTPUT - CODE                                                      **
C* OUTPUT - I1                                                        **
C* OUTPUT - I2                                                        **
C* OUTPUT - I3                                                        **
C* OUTPUT - I4                                                        **
C* OUTPUT - F1                                                        **
C* OUTPUT - F2                                                        **
C* OUTPUT - F3                                                        **
C* OUTPUT - F4                                                        **
C* OUTPUT - F5                                                        **
C* OUTPUT - F6                                                        **
C* INPUT  - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       PARSIT                                                 **
C* called by   NEC2D                                                  **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL PARSIT

C     Dummy arguments.
      CHARACTER*(*) CODE
      INTEGER INUNIT , I1 , I2 , I3 , I4 , IFAIL , DEBUG
      REAL*8 F1 , F2 , F3 , F4 , F5 , F6

C     Local variables.
      INTEGER IEOF
      INTEGER INTVAL(4)
      REAL*8  REAVAL(6)
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING READMN'
      ENDIF

C  Call the routine to read the record and parse it.

      CALL PARSIT(INUNIT,4,6,CODE,INTVAL,REAVAL,IEOF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN

C  Set the return variables to the buffer array elements.
      IF ( IEOF.LT.0 ) CODE = 'EN'
      I1 = INTVAL(1)
      I2 = INTVAL(2)
      I3 = INTVAL(3)
      I4 = INTVAL(4)
      F1 = REAVAL(1)
      F2 = REAVAL(2)
      F3 = REAVAL(3)
      F4 = REAVAL(4)
      F5 = REAVAL(5)
      F6 = REAVAL(6)
 
      RETURN
 
      END
