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

add_executable( aegnec2 cnec2d.c cnec2d.h nec.c nec.h )
target_link_libraries( aegnec2 nec2d ${NEC_LIBS} )
set_target_properties( aegnec2 PROPERTIES LINKER_LANGUAGE Fortran )
install( TARGETS aegnec2 RUNTIME DESTINATION bin )

