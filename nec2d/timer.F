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

#include <config.h>
 
C*--------------------------------------------------------------------**

      SUBROUTINE TIMER( CPUSECS )

C*--------------------------------------------------------------------**
C*                                                                    **
C* TIMER returns cpu time in seconds.  Must be customized!!!          **
C*                                                                    **
C*------------------------- Dummy Arguments --------------------------**
C*                                                                    **
C* OUTPUT - CPUSECS                                                   **
C*                                                                    **
C*------------------------- COMMON Variables -------------------------**
C*                                                                    **
C* modifies    ** NOTHING **                                          **
C* uses value  ** NOTHING **                                          **
C*                                                                    **
C*----------------------- External Subprograms -----------------------**
C*                                                                    **
C* calls       ** NOTHING **                                          **
C* called by   FACIO   NEC2                                           **
C*                                                                    **
C*--------------------------------------------------------------------**
 
      IMPLICIT NONE

C     External routines.
#if defined (HAVE_ETIME)

C     EXTERNAL  ETIME

#elif defined(HAVE_CPU_TIME)

C      EXTERNAL CPU_TIME

#elif defined(HAVE_SECOND)

C      EXTERNAL SECOND

#endif

C     Dummy arguments.
      REAL*8 CPUSECS

#if defined (HAVE_ETIME)

      REAL ETIME
      REAL T1, TARRAY(2)

      T1 = ETIME(TARRAY)
      CPUSECS=TARRAY(1)

#elif defined(HAVE_CPU_TIME)

       CALL CPU_TIME(CPUSECS)

#elif defined(HAVE_SECOND)

       CALL SECOND(CPUSECS)

#else

      CPUSECS=0.0D0

#endif

      RETURN
 
      END

