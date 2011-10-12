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
 * csomnec.c -- C driver routine for SOMNEC.
 *
 * This file contains the functions for driving SOMNEC from C.
 *
 */

#include <config.h>
#include <fcint.h>
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "csomnec.h"

/* 
 * Prototype for FORTRAN routine assuming one underscore and length of character 
 * array is passed at end of arg list by value as an unsigned long. This is not
 * portable but works on many machines with many compilers.
 */
void FC_somnec( INTEGER* , INTEGER* , INTEGER* , INTEGER*, REAL_8*, REAL_8*, 
               REAL_8*, char *, FTNLEN );


/*
 * Driver routine.
 */

int 
csomnec( int verbose , int debug , int ipt , double epr , double sig , 
         double fmhz , char *otfile )
{

  FTNLEN otfile_len;
  INTEGER VERBSE , DEBUG , IPT , IFAIL;
  REAL_8 EPR , SIG , FMHZ;

  /* Cast variables to Fortran types */
  VERBSE = (INTEGER) verbose;
  DEBUG  = (INTEGER) debug;
  IPT    = (INTEGER) ipt;
  EPR    = (REAL_8)epr;
  SIG    = (REAL_8)sig;
  FMHZ   = (REAL_8)fmhz;

  /* Length of job name needed for call to FORTRAN routine. */
  otfile_len = (FTNLEN) strlen( otfile );

  /* Call SOMNEC */
  IFAIL = 0;
  FC_somnec(&IFAIL , &VERBSE , &DEBUG , &IPT , &EPR , &SIG , &FMHZ , 
            otfile , otfile_len );

  return IFAIL;

}


/*
 * Maximium length of error message.
 */

enum {
  MAX_ERR_LEN = 255
};


/*
 * Array of error messages. Position in array corresponds to position of
 * respective error symbol in the enumeration in the header file.
 */

static char err_msg[CSOMNEC_ERR_END][MAX_ERR_LEN] = {
  "successful execution",           /* IFAIL =  0 */
  "error creating output file"      /* IFAIL =  1 */
};

/*
 * Return pointer to string containing error message for a given
 * error code.
 */

char *
csomnec_err( const int err_code )
{
  
  if( err_code < CSOMNEC_ERR_END )
    {
      return err_msg[err_code];
    }
  else 
    {
      fprintf( stderr, "*** Internal error: error code too big - please report bug! ***" );
      exit( EXIT_FAILURE );
    }
  
}
