#
# aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
# Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>
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

if( UNIX )

   # Test suite only works on Unix at moment

   set( RUN_NEC ${CNEC2_BINARY_DIR}/tests/runnec.sh )
   set( RUN_SOMNEC ${CNEC2_BINARY_DIR}/tests/runsomnec.sh )

   configure_file( ${CNEC2_SOURCE_DIR}/tests/runnec.sh.cmake
                   ${CNEC2_BINARY_DIR}/tests/runnec.sh )

   configure_file( ${CNEC2_SOURCE_DIR}/tests/runsomnec.sh.cmake
                   ${CNEC2_BINARY_DIR}/tests/runsomnec.sh )

    file( WRITE "${CNEC2_BINARY_DIR}/tests/timing.txt" "
------------------------------------------------
Test          Status    Fill    Factor    Run
                        Time     Time     Time
------------------------------------------------
" )

   add_test( test0300   ${BASH_PROGRAM} ${RUN_NEC} test0300 )
   add_test( test0600   ${BASH_PROGRAM} ${RUN_NEC} test0600 )
   add_test( test1200   ${BASH_PROGRAM} ${RUN_NEC} test1200 )
   add_test( test2000   ${BASH_PROGRAM} ${RUN_NEC} test2000 )
   add_test( testcomm   ${BASH_PROGRAM} ${RUN_NEC} testcomm )
   add_test( testcoup   ${BASH_PROGRAM} ${RUN_NEC} testcoup )
   add_test( testdip1   ${BASH_PROGRAM} ${RUN_NEC} testdip1 )
   add_test( testdip2   ${BASH_PROGRAM} ${RUN_NEC} testdip2 )
   add_test( testnetw   ${BASH_PROGRAM} ${RUN_NEC} testnetw )
   add_test( testsymm   ${BASH_PROGRAM} ${RUN_NEC} testsymm )
   add_test( testtran   ${BASH_PROGRAM} ${RUN_NEC} testtran )
   add_test( createngf1 ${BASH_PROGRAM} ${RUN_NEC} createngf testngf1 )
   add_test( createngf2 ${BASH_PROGRAM} ${RUN_NEC} createngf testngf2 )
   add_test( testngf1   ${BASH_PROGRAM} ${RUN_NEC} testngf1 )
   add_test( testngf2   ${BASH_PROGRAM} ${RUN_NEC} testngf2 )

  if( NOT BUILTIN_SOMNEC )

     add_test( csomnec ${BASH_PROGRAM} ${RUN_SOMNEC} testvdip 10.0 20.0 0.015 )

  else( NOT BUILTIN_SOMNEC )

     add_test( stestdip ${BASH_PROGRAM} ${RUN_NEC} stestdip )

  endif( NOT BUILTIN_SOMNEC )

  # Requires test csomnec to have run if somnec not built-in!
  add_test( testvdip ${BASH_PROGRAM} ${RUN_NEC} testvdip )

endif( UNIX )
