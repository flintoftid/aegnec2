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

      SUBROUTINE READGM( INUNIT , CODE , I1 , I2 , RR1 , RR2 , R3 , R4 ,
     &                   R5 , R6 , R7 , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* READGM reads a geometry record and parses it.                      **
C*                                                                    **
C*    CODE        two letter mnemonic code                            **
C*    I1 - I2     integer values from record                          **
C*    RR1 - R7     real values from record                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* PASSED - INUNIT                                                    **
C* OUTPUT - CODE                                                      **
C* OUTPUT - I1                                                        **
C* OUTPUT - I2                                                        **
C* OUTPUT - RR1                                                       **
C* OUTPUT - RR2                                                       **
C* OUTPUT - R3                                                        **
C* OUTPUT - R4                                                        **
C* OUTPUT - R5                                                        **
C* OUTPUT - R6                                                        **
C* OUTPUT - R7                                                        **
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
C* called by   DATAGN                                                 **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL PARSIT

C     Dummy arguments.
      CHARACTER*(*) CODE
      INTEGER INUNIT , I1 , I2 , IFAIL , DEBUG
      REAL*8 RR1 , RR2 , R3 , R4 , R5 , R6 , R7

C     Local variables.
      INTEGER IEOF
      INTEGER INTVAL(2)
      REAL*8 REAVAL(7)


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING READGM'
      ENDIF

C     Call the routine to read the record and parse it.

      CALL PARSIT(INUNIT,2,7,CODE,INTVAL,REAVAL,IEOF,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN

C     Set the return variables to the buffer array elements.

      IF ( IEOF.LT.0 ) CODE = 'GE'
      I1 = INTVAL(1)
      I2 = INTVAL(2)
      RR1 = REAVAL(1)
      RR2 = REAVAL(2)
      R3 = REAVAL(3)
      R4 = REAVAL(4)
      R5 = REAVAL(5)
      R6 = REAVAL(6)
      R7 = REAVAL(7)
 
      RETURN
 
      END
