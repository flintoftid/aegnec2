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

if( NOT BUILTIN_SOMNEC )

   add_executable( aegsomnec cdrive.c csomnec.c cdrive.h csomnec.h )
   target_link_libraries( aegsomnec somnec somtest )
   set_target_properties( aegsomnec PROPERTIES LINKER_LANGUAGE Fortran )

   install( TARGETS aegsomnec RUNTIME DESTINATION bin )

endif( NOT BUILTIN_SOMNEC )
