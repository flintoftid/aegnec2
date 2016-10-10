/*
 * cnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
 * Copyright (C) 1998-2016 Ian David Flintoft <ian.flintoft@googlemail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CONFIG_H
#define _CONFIG_H

/* Name of package of package */
#define PACKAGE ${CNEC2_PACKAGE}

/* Version number of package */
#define VERSION ${CNEC2_VERSION_STRING}

/* Define if builtin somnec support requested */
#cmakedefine BUILTIN_SOMNEC

/* Define if Fortran compiler supports CPU_TIME intrinsic. */
#cmakedefine HAVE_CPU_TIME

/* Define if Fortran compiler supports ETIME intrinsic. */
#cmakedefine HAVE_ETIME

/* Define if Fortran compiler supports SECOND intrinsic. */
#cmakedefine HAVE_SECOND

/* Define if compiling against BLAS routines. */
#cmakedefine USE_BLAS

/* Define if Fortran compiler supports FLUSH intrinsic. */
#cmakedefine HAVE_FLUSH

#endif /* ifndef _CONFIG_H */
