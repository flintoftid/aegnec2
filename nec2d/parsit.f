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

      SUBROUTINE PARSIT( INUNIT , MAXINT , MAXREA , CMND , INTFLD ,
     &                   REAFLD , IEOF , IFAIL , DEBUG )

C*--------------------------------------------------------------------**
C*                                                                    **
C* PARSIT reads an input record and parses it.                        **
C*                                                                    **
C* UPDATED:  21 July 87                                               **
C*                                                                    **
C* Passed variables                                                   **
C*    MAXINT     total number of integers in record                   **
C*    MAXREA     total number of real values in record                **
C*    CMND       two letter mnemonic code                             **
C*    INTFLD     integer values from record                           **
C*    REAFLD     real values from record                              **
C*                                                                    **
C* Internal Variables                                                 **
C*    BGNFLD     list of starting indices                             **
C*    BUFFER     text buffer                                          **
C*    ENDFLD     list of ending indices                               **
C*    FLDTRM     flag to indicate that pointer is in field position   **
C*    REC        input line as read                                   **
C*    TOTCOL     total number of columns in REC                       **
C*    TOTFLD     number of numeric fields                             **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* INPUT  - INUNIT                                                    **
C* INPUT  - MAXINT                                                    **
C* INPUT  - MAXREA                                                    **
C* OUTPUT - CMND                                                      **
C* OUTPUT - INTFLD                                                    **
C* OUTPUT - REAFLD                                                    **
C* OUTPUT - IEOF                                                      **
C* OUTPUT - IFAIL                                                     **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       UPCASE                                                 **
C* called by   READGM  READMN                                         **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
      EXTERNAL UPCASE

C     Parameter definitions.
      INCLUDE 'nec2d.inc'

C     Dummy arguments.
      INTEGER INUNIT , MAXINT , MAXREA , IEOF , IFAIL , DEBUG
      INTEGER INTFLD(MAXINT)
      REAL*8 REAFLD(MAXREA)
      CHARACTER CMND*2

C     Local variables.
      LOGICAL FLDTRM, POSCOM
      CHARACTER BUFFER*20 , REC*160
      CHARACTER TBUFF*20
      CHARACTER*1 COMMTS(3)
      INTEGER BGNFLD(12) , ENDFLD(12) , TOTCOL , TOTFLD
      INTEGER I , IND , INDE , J , K , LAST , LENGTH 
 
C     Data initialisation.
      DATA COMMTS /'#', '%', ';'/


C     IDF - zero fields.
      DO 1000 I = 1,MAXINT
         INTFLD(I) = 0
 1000 CONTINUE
      DO 1010 I = 1,MAXREA
         REAFLD(I) = 0.0D0
 1010 CONTINUE

      READ (INUNIT,8000,IOSTAT=IEOF) REC
      CALL UPCASE(REC,REC,TOTCOL,DEBUG)


      IF (DEBUG.GT.0) THEN
         PRINT*, 'DEBUG: ENTERING PARSIT'
      ENDIF

C  Store opcode and clear field arrays.
      POSCOM=.FALSE.
      CMND = REC(1:2)

C  IDF Check for obvious comment line.
      DO 2010 I = 1, 3
         IF(COMMTS(I).EQ.REC(1:1)) THEN
            CMND='IG'
            RETURN
         ENDIF
 2010 CONTINUE
      IF( REC(1:1).EQ.' ' ) THEN
         DO 2012 I = 1, 3
            IF(COMMTS(I).EQ.REC(2:2)) THEN
               CMND='IG'
               RETURN
            ENDIF
 2012    CONTINUE
      ENDIF
      IF(CMND.EQ.'  ') POSCOM=.TRUE.
C  end IDF.
      DO 3000 I = 1 , MAXINT
         INTFLD(I) = 0
 3000 CONTINUE
      DO 3010 I = 1 , MAXREA
         REAFLD(I) = 0.0D0
 3010 CONTINUE
      DO 3020 I = 1 , 12
         BGNFLD(I) = 0
         ENDFLD(I) = 0
 3020 CONTINUE

C  Find the beginning and ending of each field as well as the total 
C  number of fields.

      TOTFLD = 0
      FLDTRM = .FALSE.
      LAST = MAXREA + MAXINT
      DO 4000 J = 3 , TOTCOL
         K = ICHAR(REC(J:J))

C  Check for end of line comment (`!' or '#').  This is a new 
C  modification to allow VAX-like comments at the end of data records, 
C  i.e.
C       GW 1 7 0 0 0 0 0 .5 .0001 ! DIPOLE WIRE
C       GE # END OF GEOMETRY

         IF ( K.EQ.33 .OR. K.EQ.35 ) THEN
            IF ( FLDTRM ) ENDFLD(TOTFLD) = J - 1
            GOTO 5000

C  Set the ending index when the character is a comma or space and the
C  pointer is in a field position (FLDTRM = .TRUE.).

         ELSEIF ( K.EQ.32 .OR. K.EQ.44 ) THEN
            IF ( FLDTRM ) THEN
               ENDFLD(TOTFLD) = J - 1
               FLDTRM = .FALSE.
            ENDIF

C  Set the beginning index when the character is not a comma or space
C  and the pointer is not currently in a field position (FLDTRM = .
C  FALSE).

         ELSEIF ( .NOT.FLDTRM ) THEN
            TOTFLD = TOTFLD + 1
            FLDTRM = .TRUE.
            BGNFLD(TOTFLD) = J
         ENDIF
 4000 CONTINUE
      IF ( FLDTRM ) ENDFLD(TOTFLD) = TOTCOL
 
C  Check to see if the total number of value fields is within the 
C  precribed limits.
 
 5000 IF ( TOTFLD.EQ.0 ) THEN
C IDF other possible way to get comment.
         IF ( POSCOM.EQV..TRUE. ) THEN
            CMND='IG'
         ENDIF
C end IDF.
         RETURN
      ELSEIF ( TOTFLD.GT.LAST ) THEN
         WRITE (CHRSLT,8001)
         GOTO 9010
      ENDIF
      J = MIN(TOTFLD,MAXINT)
 
C  Parse out integer values and store into integer buffer array.
 
      DO 5090 I = 1 , J
         LENGTH = ENDFLD(I) - BGNFLD(I) + 1
         BUFFER = REC(BGNFLD(I):ENDFLD(I))
         IND = INDEX(BUFFER(1:LENGTH),'.')
         IF ( IND.GT.0 .AND. IND.LT.LENGTH ) GOTO 9000
         IF ( IND.EQ.LENGTH ) LENGTH = LENGTH - 1
         READ (BUFFER(1:LENGTH),*,ERR=9000) INTFLD(I)
 5090 CONTINUE
 
C  Parse out real values and store into real buffer array.
 
      IF ( TOTFLD.GT.MAXINT ) THEN
         J = MAXINT + 1
         DO 6000 I = J , TOTFLD
            LENGTH = ENDFLD(I) - BGNFLD(I) + 1
            BUFFER = REC(BGNFLD(I):ENDFLD(I))
            IND = INDEX(BUFFER(1:LENGTH),'.')
            IF ( IND.EQ.0 ) THEN
               INDE = INDEX(BUFFER(1:LENGTH),'E')
               LENGTH = LENGTH + 1
               IF ( INDE.EQ.0 ) THEN
                  BUFFER(LENGTH:LENGTH) = '.'
               ELSE
C  Kludge to fix bug in PGI compiler.
C                 BUFFER = BUFFER(1:INDE-1)//'.'//BUFFER(INDE:LENGTH-1)
                  TBUFF = BUFFER(1:INDE-1)//'.'//BUFFER(INDE:LENGTH-1)
                  BUFFER = TBUFF
C  End kludge.
               ENDIF
            ENDIF
            READ (BUFFER(1:LENGTH),*,ERR=9000) REAFLD(I-MAXINT)
 6000    CONTINUE
      ENDIF
      RETURN

C  Print out text of record line when error occurs.
 
 9000 IF ( I.LE.MAXINT ) THEN
         WRITE (CHRSLT,8002) I
      ELSE
         I = I - MAXINT
         WRITE (CHRSLT,8003) I
      ENDIF
 9010 WRITE (CHRSLT,8004) REC
      IFAIL=28
      RETURN
 
 8000 FORMAT (A160)
 8001 FORMAT (//,' ***** CARD ERROR - TOO MANY FIELDS IN RECORD')
 8002 FORMAT (//,' ***** CARD ERROR - INVALID NUMBER AT INTEGER',
     &        ' POSITION ',I1)
 8003 FORMAT (//,' ***** CARD ERROR - INVALID NUMBER AT REAL',
     &        ' POSITION ',I1)
 8004 FORMAT (' ***** TEXT -->  ',A160)
 
      END
