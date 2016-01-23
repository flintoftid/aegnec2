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

      PROGRAM FDRIVE

C*--------------------------------------------------------------------**
C*                                                                    **
C* PROGRAM TO GENERATE NEC INTERPOLATION GRIDS FOR FIELDS DUE TO      **
C* GROUND.  FIELD COMPONENTS ARE COMPUTED BY NUMERICAL EVALUATION     **
C* OF MODIFIED SOMMERFELD INTEGRALS.                                  **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       SOMNEC                                                 **
C* called by   ** NOTHING **                                          **
C*                                                                    **
C*--------------------------------------------------------------------**

      IMPLICIT NONE

C     External routines.
      EXTERNAL SOMNEC

C     Local variables.
      CHARACTER*32 OTFILE
      INTEGER IFAIL , VERBSE , DEBUG , IPT 
      REAL*8 EPR , FMHZ , SIG


C     Read ground parameters - EPR = relative dielectric constant
C                              SIG = conductivity (mhos/m)
C                              FMHZ = frequency (MHz)
C                              IPT = 1 to print grids, = 0 otherwise.
C     IF SIG .LT. 0. then complex dielectric constant = EPR + J*SIG
C     and FMHZ is not used.

C     READ 15, EPR,SIG,FMHZ,IPT
      PRINT 100
      PRINT 101
      PRINT 102
      PRINT 103
      READ *, EPR
      PRINT 104
      READ *, SIG
      PRINT 105
      READ *, FMHZ
      PRINT 106
      READ *, IPT
      PRINT 107
      READ 24, OTFILE
      PRINT *, ' Relative Dielectric Constant = ', EPR
      PRINT *, ' Conductivity (Mhos/Meter) = ', SIG
      PRINT *, ' Frequency, MHz = ', FMHZ
      PRINT *, ' Printing Flag = ', IPT
      PRINT *, ' Data Output File Name = ', OTFILE

      VERBSE = 1
      DEBUG = 0
      IFAIL = 0

      CALL SOMNEC( IFAIL , VERBSE , DEBUG , IPT , EPR , SIG , FMHZ , 
     &             OTFILE )

 100  FORMAT (' Program to Calculate Ground Interpolation Grid')
 101  FORMAT (' For NEC2 Using Sommerfeld-Norton Method')
 102  FORMAT (' ')
 103  FORMAT (' Enter Relative Dielectric Constant:')
 104  FORMAT (' Enter Conductivity (Mhos/Meter):')
 105  FORMAT (' Enter Frequency (MHz):')
 106  FORMAT (' Enter 1 to Print Grids, 0 to Suppress Printing:')
 107  FORMAT (' Enter Data Output Filename:')
 24   FORMAT (A)

      END
