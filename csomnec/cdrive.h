/*
 * cnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
 * Copyright (C) 1998-2009 Ian David Flintoft <idf1@ohm.york.ac.uk>
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
 * cdrive.h -- Header file for C front end to SOMNEC.
 *
 * Contains the arrays bounds and other parameters.
 *
 */

/* 
 * NEC input file suffix. This must be same as SUFINP in nec2d.inc 
 */

#define SUFSOM          "som"

/* 
 * Default sizes and bounds 
 */

enum {

  DEBUG_MIN     = 0,
  DEBUG_DEFAULT = 0,
  DEBUG_MAX     = 10,         

};

/* 
 * Prototypes 
 */

void get_options( int, char **, char **, int *, int *, int *, int *, 
		  double *, double *, double * );

void get_pos_long( const char, long *, char *, const long, const long );

void get_pos_float( const char, double *, char * );

void version_info( int );

void usage( const int, const char*, const char* );

void verify_output_file( char ** , char ** );

