#!/bin/sh
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
# ---------------------------------------------------------------------
# 
# runsomnec.sh(.cmake) -- Control script for cnec2 testsuite.
# 
# ---------------------------------------------------------------------

# Script name.
script=`basename $0`

# Top of source tree.    
top_srcdir="@CNEC2_SOURCE_DIR@"

# Top of build tree.    
top_bindir="@CNEC2_BINARY_DIR@"

# The nec executable.
CSOMNEC="@CNEC2_BINARY_DIR@/csomnec/aegsomnec"

# Test source dir
test_srcdir="@CNEC2_SOURCE_DIR@/tests"

# Test run dir
test_bindir="@CNEC2_BINARY_DIR@/tests"

# Return code from script.
ret=0

# How to suppres newline in echo.
ECHO_N=""
ECHO_C=""
ECHO_T=""

# Check we have a valid cnec2 executable.
if test ! -x "$CSOMNEC"
then
   echo "$script: cannot execute $CSOMNEC"   
   exit 1
fi
if ( ( eval $CSOMNEC --version ) 2>&1 >/dev/null )
then
   :
else
   echo "$script: invalid csomnec executable: $CSOMNEC"   
   exit 1
fi

testname="$1"
freq="$2"
epsr="$3"
sigma="$4"

echo "Test file $testname"
echo

# Run the test, careful may be VPATH build.
echo $ECHO_N "    running test... $ECHO_C"
"$CSOMNEC" -f "$freq" -e "$epsr" -c "$sigma" "$testname.som" 2>&1 > "$testname.log"
echo "$ECHO_T done."

echo

exit "$ret"
