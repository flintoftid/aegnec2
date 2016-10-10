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
 
      SUBROUTINE PRNT( IN1 , IN2 , IN3 , FL1 , FL2 , FL3 , FL4 , FL5 ,
     &                 FL6 , CTYPE , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C*  PRNT prints the input data for impedance loading, inserting blanks**
C*  for numbers that are zero.                                        **
C*                                                                    **
C*  INPUT:                                                            **
C*  IN1-3 = INTEGER VALUES TO BE PRINTED                              **
C*  FL1-6 = REAL VALUES TO BE PRINTED                                 **
C*  CTYPE = CHARACTER STRING TO BE PRINTED                            **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - IN1                                                       **
C* INPUT  - IN2                                                       **
C* INPUT  - IN3                                                       **
C* INPUT  - FL1                                                       **
C* INPUT  - FL2                                                       **
C* INPUT  - FL3                                                       **
C* INPUT  - FL4                                                       **
C* INPUT  - FL5                                                       **
C* INPUT  - FL6                                                       **
C* INPUT  - CTYPE                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   LOAD                                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Duumy arguments.
      INTEGER IN1 , IN2 , IN3 , DEBUG
      REAL*8 FL1 , FL2 , FL3 , FL4 , FL5 , FL6
      CHARACTER CTYPE*(*)

C     Local variables.
      INTEGER I
      CHARACTER CINT(3)*5 , CFLT(6)*13
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING PRNT'
      ENDIF
 
      DO 1 I = 1 , 3
         CINT(I) = '     '
    1 CONTINUE
      IF ( IN1.EQ.0 .AND. IN2.EQ.0 .AND. IN3.EQ.0 ) THEN
         CINT(1) = '  ALL'
      ELSE
         IF ( IN1.NE.0 ) WRITE (CINT(1),90) IN1
         IF ( IN2.NE.0 ) WRITE (CINT(2),90) IN2
         IF ( IN3.NE.0 ) WRITE (CINT(3),90) IN3
      ENDIF
      DO 2 I = 1 , 6
         CFLT(I) = '     '
    2 CONTINUE
      IF ( DABS(FL1).GT.1.0D-30 ) WRITE (CFLT(1),91) FL1
      IF ( DABS(FL2).GT.1.0D-30 ) WRITE (CFLT(2),91) FL2
      IF ( DABS(FL3).GT.1.0D-30 ) WRITE (CFLT(3),91) FL3
      IF ( DABS(FL4).GT.1.0D-30 ) WRITE (CFLT(4),91) FL4
      IF ( DABS(FL5).GT.1.0D-30 ) WRITE (CFLT(5),91) FL5
      IF ( DABS(FL6).GT.1.0D-30 ) WRITE (CFLT(6),91) FL6
      WRITE (CHRSLT,92) (CINT(I),I=1,3) , (CFLT(I),I=1,6) , CTYPE
 
      RETURN
 
   90 FORMAT (I5)
   91 FORMAT (1P,E13.4)
   92 FORMAT (/,3X,3A,3X,6A,3X,A)
 
      END
