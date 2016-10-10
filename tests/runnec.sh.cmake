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
# runnec.sh(.cmake) -- Control script for cnec2 testsuite.
# 
# ---------------------------------------------------------------------


# Script name.
script=`basename $0`

# Top of source tree.    
top_srcdir="@CNEC2_SOURCE_DIR@"

# Top of build tree.    
top_bindir="@CNEC2_BINARY_DIR@"

# The nec executable.
CNEC2="@CNEC2_BINARY_DIR@/cnec2/aegnec2"

# Name of results timing file.
timing_file="@CNEC2_BINARY_DIR@/tests/timing.txt"

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
if test ! -x "$CNEC2"
then
   echo "$script: cannot execute $CNEC2"   
   exit 1
fi
if ( ( eval $CNEC2 --version ) 2>&1 >/dev/null )
then
   :
else
   echo "$script: invalid cnec2 executable: $CNEC2"   
   exit 1
fi

testname="$1"
copyname="$2"

valid_file="@CNEC2_SOURCE_DIR@/tests/$testname.val"

test_file="$testname.nec"
result_file="$testname.res"
result_file_comp="$testname.res.comp"
valid_file_comp="$testname.val.comp"
log_file="$testname.log"

echo "Test file $testname"
echo

# Run the test, careful may be VPATH build.
/bin/cp -f "@CNEC2_SOURCE_DIR@/tests/$test_file" .
/bin/rm -f "$result_file"

echo $ECHO_N "    running test... $ECHO_C"
"$CNEC2" -s 2000 "$testname" 2>&1 > "$log_file"
echo "$ECHO_T done."

# Determine fill, factor and run times. Average if multiple
# frequency points are used.

temp=`cat $result_file    | grep FACTOR           \
                          | sed -e "s/  */ /g"    \
                          | cut -d" " -f3`

if test -z "$temp"
then
   fill=0.0
   factor=0.0
else
   fill=`cat $result_file    | grep FACTOR           \
                             | sed -e "s/  */ /g"    \
		             | cut -d" " -f3         \
		             | nawk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/NR}'`

   factor=`cat $result_file  | grep FACTOR           \
                             | sed -e "s/  */ /g"    \
                             | cut -d" " -f6         \
                             | nawk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/NR}'`
fi

runtime=`cat $result_file | grep RUN              \
                          | sed -e "s/  */ /g"    \
                          | cut -d" " -f5`

echo "    fill time = $fill, factor time = $factor, run time = $runtime"

# If validation file exists for test check results against it.
if test -f "$valid_file"
then
        
   echo $ECHO_N "    validating results... $ECHO_C"
   cat "$result_file" | grep -v "RUN TIME" | grep -v "FACTOR=" \
                      | sed -e "s/ -0.00 /  0.00 /g" \
                      | sed -e "s/ -0.0000 /  0.0000 /g" \
                      | sed -e "s/ -0.00000 /  0.00000 /g" \
                      > "$result_file_comp"
   cat "$valid_file"  | grep -v "RUN TIME" | grep -v "FACTOR="  \
                      | sed -e "s/ -0.00 /  0.00 /g" \
                      | sed -e "s/ -0.0000 /  0.0000 /g" \
                      | sed -e "s/ -0.00000 /  0.00000 /g" \
                      > "$valid_file_comp"

   diffs=`diff $result_file_comp $valid_file_comp`

   if test -z "$diffs"
   then
      echo "$ECHO_T ok."
      pass="pass"
   else
      echo "$ECHO_T failed. See $log_file."
      diff "$result_file_comp" "$valid_file_comp" >> "$log_file"
      pass="fail"
      ret=1
   fi

else

   echo "    no validation file."
   pass="????"

fi

echo "$testname  $pass  $fill  $factor  $runtime" | nawk '{printf("%-14s %4s %8.3f %8.3f %8.3f\n",$1,$2,$3,$4,$5)}' >> "$timing_file"

rm -f "$result_file_comp"
rm -f "$valid_file_comp"

if [ -n "$copyname" ]
then

   for suffix in res log imp ngf som 
   do

      /bin/mv -f "$testname.$suffix" "$copyname.$suffix"
       echo "$testname.$suffix" "$copyname.$suffix"

   done

fi

echo

exit "$ret"
