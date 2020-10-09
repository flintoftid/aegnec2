#
# aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
# Copyright (C) 1998-2016 Ian David Flintoft <ian.flintoft@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

set( CNEC2_MANPAGES "aegnec2.1" )

if( NOT BUILTIN_SOMNEC )
   set( CNEC2_MANPAGES ${CNEC2_MANPAGES} "aegsomnec.1" )
endif( NOT BUILTIN_SOMNEC )
 
file( WRITE ${CNEC2_BINARY_DIR}/doc/License.txt "
aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2
Copyright (C) 1998-2016 Ian David Flintoft <ian.flintoft@googlemail.com>

This version of NEC2 is distributed under the terms of the GNU General
Public License Version 3. See below for full license.
" )

if( USE_MLAPACK )
   file( APPEND ${CNEC2_BINARY_DIR}/doc/License.txt "
This build of cnec includes software from the LAPACK 3.2.1
mathematics libraries. These files are distributed under their own GPL
Version 3 compatible license. See below for full license.
" )
endif( USE_MLAPACK )

file( APPEND ${CNEC2_BINARY_DIR}/doc/Licence.txt "\n\nGNU General Public License Version 3\n------------------------------------\n\n" )
file( READ "gpl-3.0.txt" TEXT_VAR )
file( APPEND ${CNEC2_BINARY_DIR}/Licence.txt ${TEXT_VAR} )
if( USE_MLAPACK )
   file( APPEND ${CNEC2_BINARY_DIR}/Licence.txt "\n\nLAPACK License\n--------------\n\n" )
   file( READ "${CNEC2_SOURCE_DIR}/mlapack/LAPACK_LICENSE.txt" TEXT_VAR )
   file( APPEND ${CNEC2_BINARY_DIR}/doc/Licence.txt ${TEXT_VAR} )
endif( USE_MLAPACK )

install( FILES ${CNEC2_MANPAGES} DESTINATION man/man1 )

install( FILES Bugs.md ${CNEC2_SOURCE_DIR}/ChangeLog.txt 
         ${CNEC2_BINARY_DIR}/doc/License.txt News.md ToDo.md DESTINATION doc/aegnec2 )

