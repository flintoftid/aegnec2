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
 
      SUBROUTINE TRIO( J , LD , ICON1 , ICON2 , SI , BI , JMAX , 
     &                 JSNO , JCO , AX , BX , CX , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* Compute the components of all basis functions on segment J due to  **
C* each of the segments connected to J.                               **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - J                                                         **
C* INPUT  - LD                                                        **
C* INPUT  - ICON1                                                     **
C* INPUT  - ICON2                                                     **
C* PASSED - SI                                                        **
C* PASSED - BI                                                        **
C* INPUT  - JMAX                                                      **
C* OUTPUT - JSNO                                                      **
C* OUTPUT - JCO                                                       **
C* PASSED - AX                                                        **
C* PASSED - BX                                                        **
C* PASSED - CX                                                        **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       SBF                                                    **
C* called by   CMNGF   CMSET   CMSW                                   **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL SBF

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER J , LD , JMAX , JSNO , IFAIL , DEBUG
      INTEGER ICON1(2*LD) , ICON2(2*LD) , JCO(JMAX)
      REAL*8 SI(LD) , BI(LD) , AX(JMAX) , BX(JMAX) , CX(JMAX) 

C     Local variables.
      INTEGER IEND , JCOX , JEND
 

      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING TRIO'
      ENDIF
 
      JSNO = 0
      JCOX = ICON1(J)
      IF ( JCOX.GT.PATOFF ) GOTO 7
      JEND = -1
      IEND = -1
      IF ( JCOX ) 1 , 7 , 2
    1 JCOX = -JCOX
      GOTO 3
    2 JEND = -JEND
    3 IF ( JCOX.EQ.J ) GOTO 6
      JSNO = JSNO + 1
      IF ( JSNO.GE.JMAX ) GOTO 9
      CALL SBF(JCOX,J,AX(JSNO),BX(JSNO),CX(JSNO),LD,ICON1,ICON2,SI,
     &         BI,JMAX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      JCO(JSNO) = JCOX
      IF ( JEND.EQ.1 ) GOTO 4
      JCOX = ICON1(JCOX)
      GOTO 5
    4 JCOX = ICON2(JCOX)
    5 IF ( JCOX ) 1 , 9 , 2
    6 IF ( IEND.EQ.1 ) GOTO 8
    7 JCOX = ICON2(J)
      IF ( JCOX.GT.PATOFF ) GOTO 8
      JEND = 1
      IEND = 1
      IF ( JCOX ) 1 , 8 , 2
    8 JSNO = JSNO + 1
      CALL SBF(J,J,AX(JSNO),BX(JSNO),CX(JSNO),LD,ICON1,ICON2,SI,BI,
     &         JMAX,IFAIL,DEBUG)
      IF(IFAIL.NE.0) RETURN
      JCO(JSNO) = J
      RETURN
    9 WRITE (CHRSLT,10) J
 
      IFAIL=33
      RETURN
 
   10 FORMAT (' TRIO - SEGMENT CONNENTION ERROR FOR SEGMENT',I5)
 
      END
