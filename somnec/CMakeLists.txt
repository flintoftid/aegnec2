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

add_library ( somnec STATIC bessel.f evlua.f gshank.f hankel.f 
                            lambda.f rom1.f saoa.f somnec.f  cntour.inc 
                            evlcom.inc ggrid.inc )

add_library ( somtest STATIC test.f timer.F )

if( NOT BUILTIN_SOMNEC )
   add_executable( fsomnec fdrive.f )
   target_link_libraries( fsomnec somnec somtest )
   set_target_properties( fsomnec PROPERTIES LINKER_LANGUAGE Fortran )
endif( NOT BUILTIN_SOMNEC )

