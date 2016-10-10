#
# cnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
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

project ( CNEC2 C Fortran )

cmake_minimum_required( VERSION 2.8.0 )

set( CNEC2_PACKAGE "\"cnec2\"" )
set( CNEC2_VERSION_MAJOR "0" )
set( CNEC2_VERSION_MINOR "9" )
set( CNEC2_VERSION_PATCH "0" )
set( CNEC2_VERSION_STRING "\"${CNEC2_VERSION_MAJOR}.${CNEC2_VERSION_MINOR}.${CNEC2_VERSION_PATCH}\"" )

# Look here for includes and libraries.
include_directories ( ${CNEC2_BINARY_DIR} ) 
link_directories ( ${CNEC2_BINARY_DIR}/nec2d ) 

# Set up C to Fortran interface.
include( FortranCInterface )
FortranCInterface_HEADER( "${CNEC2_BINARY_DIR}/fcint.h"
                          MACRO_NAMESPACE "FC_"
                          SYMBOL_NAMESPACE "FC_"
                          SYMBOLS nec2d somnec )

FortranCInterface_VERIFY()

# Determine if to use builtin somnec functionality.
option( BUILTIN_SOMNEC  "Build somnec functionality into cnec2 programme" ON )

# Find supported Fortran timing routine.
try_compile( HAVE_ETIME ${CNEC2_BINARY_DIR}/CMakeTests ${CNEC2_SOURCE_DIR}/CMakeTests/testetime.f )
if( NOT HAVE_ETIME )
   try_compile( HAVE_CPU_TIME ${CNEC2_BINARY_DIR}/CMakeTests ${CNEC2_SOURCE_DIR}/CMakeTests/testcputime.f )
   if( NOT HAVE_CPU_TIME )
      try_compile( HAVE_SECOND ${CNEC2_BINARY_DIR}/CMakeTests ${CNEC2_SOURCE_DIR}/CMakeTests/testsecond.f )
   endif( NOT HAVE_CPU_TIME )   
endif( NOT HAVE_ETIME )

# See if Fortran compiler supports FLUSH intrinsic.
set( HAVE_FLUSH ${CNEC2_BINARY_DIR}/CMakeTests ${CNEC2_SOURCE_DIR}/CMakeTests/testflush.f )

# Determine which solver to use.
option( FIND_BLAS  "Attempt to find useable BLAS library" ON )
option( FIND_LAPACK "Attempt to find usable LAPAKC library" ON )

if( NOT SOLVER )
   message( STATUS "Finding optimum solver ..." )
   if( FIND_LAPACK )
      find_package( LAPACK )
   endif( FIND_LAPACK )
   if( LAPACK_FOUND )
     set( USE_BLAS TRUE )
     set( USE_MLAPACK FALSE )
   else( LAPACK_FOUND )
      if( FIND_BLAS )
         find_package( BLAS )
      endif( FIND_BLAS )
      if( BLAS_FOUND )
         set( USE_BLAS TRUE )
         set( USE_MLAPACK TRUE )
      endif( BLAS_FOUND )
   endif( LAPACK_FOUND)
   message( STATUS "Building for solver: ${SOLVER}" )
elseif( SOLVER STREQUAL "BUILTIN" )
   message( STATUS "Building for solver: ${SOLVER}" )
   set( USE_BLAS FALSE )
   set( USE_MLAPACK FALSE )
elseif( SOLVER STREQUAL "BLAS" )
   message( STATUS "Building for solver: ${SOLVER}" )
   if( FIND_BLAS )
         find_package( BLAS REQUIRED )
      endif( FIND_BLAS )
   set( USE_BLAS TRUE )
   set( USE_MLAPACK TRUE )
elseif( SOLVER STREQUAL "LAPACK" )
   message( STATUS "Building for solver: ${SOLVER}" )
   if( FIND_LAPACK )
      find_package( LAPACK REQUIRED )
   endif( FIND_LAPACK )
   set( USE_BLAS TRUE )
   set( USE_MLAPACK FALSE )
else( NOT SOLVER)
   message( FATAL_ERROR "Invalid solver: " SOLVER ) 
endif( NOT SOLVER )

# Set nec link libraries according to solver and somnec options.
set( NEC_LIBS )
if( BUILTIN_SOMNEC )
   set( NEC_LIBS ${NEC_LIBS} somnec )
endif( BUILTIN_SOMNEC )
if( USE_BLAS AND USE_MLAPACK )
   set( NEC_LIBS ${NEC_LIBS} ${BLAS_LIBRARIES} mlapack )
elseif( USE_BLAS AND NOT USE_MLAPACK )
   set( NEC_LIBS ${NEC_LIBS} ${LAPACK_LIBRARIES} )
endif( USE_BLAS AND USE_MLAPACK )

# Configure global include file.
configure_file( ${CNEC2_SOURCE_DIR}/config.h.cmake
                ${CNEC2_BINARY_DIR}/config.h )

# Sub-directories.
if( USE_BLAS AND USE_MLAPACK )
   add_subdirectory ( mlapack )
endif( USE_BLAS AND USE_MLAPACK )

if( NOT BUILTIN_SOMNEC )
   add_subdirectory ( csomnec )
endif( NOT BUILTIN_SOMNEC )

add_subdirectory( nec2d )
add_subdirectory( somnec )
add_subdirectory( cnec2 )

add_subdirectory( doc )

# Package testing.
enable_testing()
find_program( BASH_PROGRAM bash )
add_subdirectory ( tests )

# Installation.
install( FILES ReadMe.txt Licence.txt DESTINATION doc/aegnec2 )

# Package generation.
include( InstallRequiredSystemLibraries )
set( CPACK_PACKAGE_NAME "aegnec2" )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "AEGNEC2 CEM code" )
set( CPACK_PACKAGE_VENDOR "Ian Flintoft <ian.flintoft@googlemail.com>" )
set( CPACK_PACKAGE_CONTACT "Ian Flintoft <ian.flintoft@googlemail.com>" )
set( CPACK_PACKAGE_DESCRIPTION_FILE "${CNEC2_SOURCE_DIR}/ReadMe.txt" )
set( CPACK_RESOURCE_FILE_LICENSE "${CNEC2_SOURCE_DIR}/Licence.txt" )
set( CPACK_PACKAGE_VERSION_MAJOR "${CNEC2_VERSION_MAJOR}" )
set( CPACK_PACKAGE_VERSION_MINOR "${CNEC2_VERSION_MINOR}" )
set( CPACK_PACKAGE_VERSION_PATCH "${CNEC2_VERSION_PATCH}" )
set( CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${CPACK_VERSION_MAJOR}.${CPACK_VERSION_MINOR}.${CPACK_VERSION_PATCH}" )
set( CPACK_SOURCE_IGNORE_FILES "/\\\\.hg/;\\\\.swp$;\\\\.#;/#;\\\\.*~" )

if( UNIX )
   set( CPACK_GENERATOR "STGZ;TGZ" )
   set( CPACK_SOURCE_GENERATOR "TGZ" )
elseif( WIN32 )
   set( CPACK_GENERATOR "NSIS" )
   set( CPACK_SOURCE_GENERATOR "ZIP" )
endif( UNIX)

include( CPack )
