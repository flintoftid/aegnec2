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
 * somnec.c -- C front end for SOMNEC.
 *
 * This file contains the main program and options parsing routines of
 * the C front end.
 *
 */

#include <config.h>

 
#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <ctype.h>

#include <error.h>
#include <getopt.h>

#include "cdrive.h"
#include "csomnec.h"

/* 
 * Program call name. 
 */

char *program_name;


/*
 * Call SOMNEC using parameters determined by command line options. 
 */

int main( int argc, char **argv )
{

  FILE *fpout;                       /* Output file pointer */
  char *job_name    = 0;             /* Job name */
  char *somnec_file = 0;             /* Output file name */
  int ret           = 0;             /* Return code from somnec */
  int verbose       = 0;             /* Verbose output if true */
  int debug         = DEBUG_DEFAULT; /* Debug level */
  int remove_old    = 0;             /* Remove old files if true */
  int print_flag    = 0;             /* Print flag */
  double freq;                       /* Frequency (MHz) */
  double epsilon_r;                  /* Relative permittivity for ground */         
  double sigma;                      /* Conductivity of ground (S/m) */         

  /* Set program name for error messages */
  program_name = argv[0];

  /* Get options and file name */
  get_options( argc, argv, &job_name, &verbose, &debug, &remove_old, &print_flag, 
               &freq, &epsilon_r, &sigma );

  /* Verify output file name */
  verify_output_file( &job_name , &somnec_file );

  /* Give info if verbose or debug */
  if( debug > 0 || verbose )
    {
      fprintf( stderr , "Frequency    = %12.6e (MHz)\n" , freq );
      fprintf( stderr , "Permittivity = %12.6e \n" , epsilon_r );
      fprintf( stderr , "Conductivity = %12.6e (S/m)\n" , sigma );
      fprintf( stderr , "Print Flag   = %12d\n" , print_flag );
      fprintf( stderr , "Remove Flag  = %12d\n" , remove_old );
      fprintf( stderr , "Output file  = %s\n" , somnec_file );
    }

  /* Remove old output files */
  if( remove_old )
    {
      fpout = fopen( somnec_file , "r" );
      if ( fpout )
	{
	  fclose( fpout );
	  ret = remove( somnec_file );
	}
    }

  ret = csomnec( verbose , debug , print_flag , epsilon_r , sigma , freq , somnec_file );
  
  if( ret != CSOMNEC_ERR_OK ) 
    {
      free( somnec_file );
      error( 1, 0, "%s", csomnec_err( ret ) );
    }
  
  free( somnec_file );
  
  return 0;

}


/*
 * Give version info and exit.
 */

void 
version_info( int verbose )
{
  
  printf( "\nsomnec (%s) %s\n", PACKAGE, VERSION );
  printf( "Copyright (C) 1999-2016 Ian David Flintoft\n" );
  printf( "Maintained by Ian D. Flintoft <ian.flintoft@googlemail.com>\n" );
  printf( "%s comes with NO WARRANTY,\n", PACKAGE );
  printf( "to the extent permitted by law.\n" );
  printf( "You may redistribute copies of %s\n", PACKAGE );
  printf( "under the terms of the GNU General Public License.\n" );
  printf( "For more information about these matters,\n" );
  printf( "see the file named COPYING.\n\n" );

  printf( "Built on %s at %s.\n\n", __DATE__, __TIME__ );

  exit( EXIT_SUCCESS );

  return;

}


/* 
 * Give usage info and exit. 
 */

void 
usage( const int exitcode, const char *error, const char *add )
{
  
  printf( "\nsomnec (%s) %s\n", PACKAGE, VERSION );
  printf( "\n" );
  printf( "Usage:  somnec [options] jobname[.som]\n\n" );
  printf( "where jobname.som is the name of a output file. The \".som\" suffix\n" );
  printf( "may be omitted. The valid options are:\n\n" );
  printf( "    -h, --help                Give this help and exit\n" );
  printf( "    -V, --version             Give version info and exit\n" );
  printf( "    -v, --verbose             Give verbose output\n" );
  printf( "    -d, --debug=N             Debug level (default 0)\n" );
  printf( "    -e, --permittivity=FLOAT  Relative permittivity of ground\n" );
  printf( "    -c, --conductivity=FLOAT  Conductivity of ground (S/m)\n" );
  printf( "    -f, --frequency=FLOAT     Frequency (MHz)\n" );
  printf( "    -r, --remove-old          Remove old output file\n" );
  printf( "    -p, --print               Print detailed table\n" );
  printf( "\nReport bugs to Ian Flintoft <ian.flintoft@googlemail.com>\n" ); 
  printf( "\n" );

  if(error) 
    fprintf( stderr, "%s: %s: %s\n", program_name, error, add );

  exit(exitcode);

  return;

}


/* 
 * Parse command line options, check bounds and determine
 * input file name.
 */

void 
get_options( int argc, char **argv , char **job_name, int *verbose, 
             int *debug, int* remove_old, int *print_flag, double *freq, 
             double *epsilon_r, double *sigma )
{
  
  static struct option long_options[] = {
    { "version",           no_argument,       0, 'V' },
    { "help",              no_argument,       0, 'h' },
    { "verbose",           no_argument,       0, 'v' },
    { "debug",             required_argument, 0, 'd' },
    { "print",             no_argument,       0, 'p' },
    { "frequency",         required_argument, 0, 'f' },
    { "conductivity",      required_argument, 0, 'c' },
    { "permittivity",      required_argument, 0, 'e' },
    { "remove-old",        no_argument,       0, 'r' },
    { 0, 0, 0, 0 }
  };
  char c;           /* Current short option */
  int got_freq = 0; /* True when frequency option found */
  int got_perm = 0; /* True when permittivity option found */
  int got_cond = 0; /* True when conductivity option found */
  long tmp_long;    /* Temporary long */

  while(1) 
    {

      c = getopt_long (argc, argv, "vhVd:pf:c:e:r", long_options, 0 );
      
      if( c == -1 ) break;
      
      switch( c ) 
	{
	case 'V':
	  version_info( *verbose );
	  exit( EXIT_SUCCESS );
	  break;
	case 'h':
	  usage( EXIT_SUCCESS, 0, 0 );
	  break;
	case 'v':
	  *verbose = 1;
	  break;
	case 'd':
	  get_pos_long( c, &tmp_long, optarg, DEBUG_MIN, DEBUG_MAX );
	  *debug = tmp_long;
	  if( *debug > 0 ) 
	    *verbose = 1;
	  break;
	case 'f':
	  get_pos_float( c, freq, optarg );
	  if( *freq <= 0.0 )
	    error( 1, 0, "value of option %c must be greater than zero", c );
	  got_freq = 1;
	  break;
	case 'c':
	  get_pos_float( c, sigma, optarg );
	  got_cond = 1;
	  break;
	case 'e':
	  get_pos_float( c, epsilon_r, optarg );
	  if( *epsilon_r < 1.0 )
	    error( 1, 0, "value of option %c must be greater than or equal to 1", c );
	  got_perm = 1;
	  break;
	case 'p':
	  *print_flag = 1;
	  break;
	case 'r':
	  *remove_old = 1;
	  break;
	case '?':
	  error( 1, 0, "ambiguous match or an extraneous parameter" );
	  exit( EXIT_FAILURE );
	  break;
	case ':':
	  error( 1, 0, "missing parameter for option" );
	  break;
	default:
	  error( 1, 0, "unknown option %c", c );
      break;
	}
      
    }
  
  if( !got_freq )
    error( 1, 0, "frequency not specified - use -f option!" );
  if( !got_perm )
    error( 1, 0, "permittivity not specified - use -e option!" );
  if( !got_cond )
    error( 1, 0, "conductivity not specified - use -c option!" );

  /* Get output file name as last arg */
  if( optind == argc )
    usage ( EXIT_FAILURE, "no output file name specified", "" );
  else if( optind == argc - 1 )
    *job_name = argv[optind];
  else
    usage ( EXIT_FAILURE, "too many non-option arguments", "" );
  
  return;

}


/*
 * Parse a field consisting of a positive integer.
 *
 */

void
get_pos_long( const char c, long *val, char *optarg, 
	      const long min, const long max )
{

  char *start;    /* Start of integer field */
  char *end;      /* End of integer field */
  char *junk;     /* Unused pointer in strtod */
  double dval;    /* Value of field as a double */
  
  printf("%s\n", optarg );

  start = optarg;
  
  if( *start == '+' ) start++;
  
  end = start;
  
  while( isdigit( *end ) ) end++;
  
  if( *end != '\0' )
    error( 1, 0, "invalid value %s for option %c", optarg, c );
  
  dval = strtod( start, &junk );
  if( dval < (long)min  || dval > (long)max )
    error( 1, 0, "value of option %c must be in range %ld to %ld", c, min, max );
  else
    *val = (long)dval;
  
  return;

}


/*
 * Parse a field consisting of a float.
 *
 */

void
get_pos_float( const char c, double *val, char *optarg )
{

  char *junk = NULL;  /* Unused pointer in strtod */
  double dval;        /* Value of field as a double */
  
  dval = strtod( optarg, &junk );

  if( junk != NULL && *junk != '\0' )
    error( 1, 0, "invalid value %s for option %c", optarg, c );

  *val = (double)dval;
  
  return;

}


/*
 * Check output file name.
 *
 */

void verify_output_file( char **job_name , char **somnec_file )
{

  char *p;              /* Pointer to start of suffix */
  char *SOMNEC_SUFFIX;  /* Input file suffix with dot */

  SOMNEC_SUFFIX = (char *) malloc( strlen( SUFSOM ) + 2 );
  strncpy( SOMNEC_SUFFIX, ".", 2 );
  strncat( SOMNEC_SUFFIX, SUFSOM, strlen( SUFSOM ) + 1 );
  p = strstr( *job_name, SOMNEC_SUFFIX );

  if( p ) 
    {
      /* Suffix in file name */
      if( strncmp ( *job_name + strlen( *job_name ) - strlen( SOMNEC_SUFFIX ), 
                    SOMNEC_SUFFIX , strlen( SOMNEC_SUFFIX ) + 1 ) != 0 ) 
	{
	  /* But not at end so ignore */
	  *somnec_file = (char *) malloc( strlen( *job_name ) + strlen( SOMNEC_SUFFIX ) + 1 );
	  strncpy( *somnec_file, *job_name, strlen( *job_name ) + 1 );
	  strncat( *somnec_file, SOMNEC_SUFFIX, strlen( SOMNEC_SUFFIX ) + 1 );
	}
    else 
      {
	/* Strip off suffix for job name */
	*somnec_file = (char *) malloc( strlen( *job_name ) + 1 );
	strncpy( *somnec_file, *job_name, strlen( *job_name ) + 1 );
	*p = '\0';
      }
    }
  else 
    {
      /* No suffix - put one on to file name */
      *somnec_file = (char *) malloc( strlen( *job_name ) + strlen( SOMNEC_SUFFIX ) + 1 );
      strncpy( *somnec_file, *job_name, strlen( *job_name ) + 1 );
      strncat( *somnec_file, SOMNEC_SUFFIX, strlen( SOMNEC_SUFFIX ) + 1 );
    }

  free( SOMNEC_SUFFIX );
  
  return;
  
}

