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
 * nec.h -- Header file for C front end to NEC2D.
 *
 * Contains the arrays bounds and other parameters.
 *
 */

/* 
 * NEC input file suffix. This must be same as SUFINP in nec2d.inc 
 */

#define SUFINP          "nec"

/* 
 * Default sizes and bounds 
 */

enum {

  MAX_LONG     = 2147400000L,
  ONE_MEGABYTE = 1048576L,
  ONE_KILOBYTE = 1024L,
  LD_MAX       = 20000L,        

  DEBUG_MIN     = 0,
  DEBUG_DEFAULT = 0,
  DEBUG_MAX     = 10,         

  NUM_SEG_MIN = 0L,
  NUM_SEG_DEFAULT = 300L,
  NUM_SEG_MAX = 10000L,         

  NUM_PATCH_MIN = 0L,
  NUM_PATCH_DEFAULT = 0L,
  NUM_PATCH_MAX = 10000L,       

  CORE_SIZE_MIN     = 0L,
  CORE_SIZE_DEFAULT = 0L,
  CORE_SIZE_MAX     = 3500L,    

  NSMAX_MIN     = 0L,
  NSMAX_DEFAULT = 30L,
  NSMAX_MAX     = 3000L,

  LOADMX_MIN     = 0L,
  LOADMX_DEFAULT = 30L,
  LOADMX_MAX     = 3000L,

  NETMX_MIN     = 0L,
  NETMX_DEFAULT = 30L,
  NETMX_MAX     = 3000L,

  FRQMAX_MIN     = 0L,
  FRQMAX_DEFAULT = 200L,
  FRQMAX_MAX     = 3000L,

  ANGMAX_MIN     = 0L,
  ANGMAX_DEFAULT = 200L,
  ANGMAX_MAX     = 3000L,

  MXCOUP_MIN     = 5L,
  MXCOUP_DEFAULT = 5L,
  MXCOUP_MAX     = 50L,

  JMAX_MIN     = 30L,
  JMAX_DEFAULT = 30L,
  JMAX_MAX     = 300L,

  NSMAXX_MIN     = 0L,
  NSMAXX_DEFAULT = 50L,
  NSMAXX_MAX     = 500L,

  NPMAX_MIN     = 0L,
  NPMAX_DEFAULT = 10L,
  NPMAX_MAX     = 500L,

  MULTI_FILE_MIN     = 0,
  MULTI_FILE_DEFAULT = 0,
  MULTI_FILE_MAX     = 2

};

/* 
 * Prototypes 
 */

void get_options( int, char**, char**, int*, double*, long*, long*, 
                  long*, long*, long*, long*, long*, long*, long*, 
		  long*, long*, int*, int* , int* );

void get_pos_long( const char, long *, char *, const long, const long );

void get_pos_bkm( const char, double *, const char *, const long, const long );

void version_info( int );

void usage( const int, const char*, const char* );

void verify_input_file( char** );

void dump_values( const char *, const int, const double, const long, 
		  const long, const long, const long, const long, 
		  const long, const long, const long, const long, 
		  const long, const long, const long, const long, 
		  const long );

void remove_old_files( const char* );




