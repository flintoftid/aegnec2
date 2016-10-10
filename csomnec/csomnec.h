/*
 * aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
 * Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>
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
 *
 * ---------------------------------------------------------------------
 *
 * csomnec.h -- C driver interface for SOMNEC.
 *
 * This file contains the interface definition for the SOMNEC C driver.
 *
 */

/* Mapping C and Fortran types. */
typedef int INTEGER;
typedef double REAL_8;
typedef unsigned long FTNLEN;

/* 
 * Return codes. Corresponding messages are in the static array in
 * csomnec.c in the same order.
 */

enum {
  CSOMNEC_ERR_OK = 0,                 /* 0 Successful execution */
  CSOMNEC_ERR_OUTFILE = 1,            /* 1 Error creating output file */ 
  CSOMNEC_ERR_END                     /* End marker */
};

/* 
 * External interface.
 */

int csomnec( int , int , int , double , double , double , char * );
char *csomnec_err( const int );

