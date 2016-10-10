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
 * nec.c -- C front end for NEC2D.
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

#include "cnec2d.h"
#include "nec.h"

/* 
 * Program call name. 
 */

char *program_name;


/*
 * Call NEC2D using the cnec2d interface with array
 * parameters determined by command line options. 
 */

int main( int argc, char **argv )
{

  char *nec_file = 0;                   /* Input file name */
  char *job_name = 0;                   /* Job name */
  int ret        = 0;                   /* Return code from cnec2d */
  int multi_file = MULTI_FILE_DEFAULT;  /* Multi file output flag */
  int verbose    = 0;                   /* Verbose output if true */
  int debug      = DEBUG_DEFAULT;       /* Debug level */
  int remove_old = 0;                   /* Remove old files if true */
  long num_seg   = NUM_SEG_DEFAULT;     /* Number of segments  */ 
  long num_patch = NUM_PATCH_DEFAULT;   /* Number of patches */
  long nsmax     = NSMAX_DEFAULT;       /* NSMAX */
  long loadmx    = LOADMX_DEFAULT;      /* LOADMX */
  long netmx     = NETMX_DEFAULT;       /* NETMX */
  long frqmax    = FRQMAX_DEFAULT;      /* Number of frequencies */
  long angmax    = ANGMAX_DEFAULT;      /* Number of angles */
  long mxcoup    = MXCOUP_DEFAULT;      /* MXCOUP */
  long jmax      = JMAX_DEFAULT;        /* JMAX */
  long nsmaxx    = NSMAXX_DEFAULT;      /* NSMAXX */
  long npmax     = NPMAX_DEFAULT;       /* NPMAX */
  long ld;                              /* LD */
  long iresrv;                          /* Core storage size multiple of 16 bytes */
  long normf;                           /* NORMF */
  double core_size = CORE_SIZE_DEFAULT; /* Core size in Mb */
  double d_num_seg;                     /* num_seg as a double */         
  double d_num_patch;                   /* num_seg as a double */         
  double d_ld;                          /* ld as a double */         
  double d_iresrv;                      /* Core storage size as a double */
  double in_core_iresrv;                /* IRESRV for in core solution */

  /* Set program name for error messages */
  program_name = argv[0];

  /* Get options and file name */
  get_options( argc, argv, &nec_file, &multi_file, &core_size, &num_seg, 
	       &num_patch, &nsmax, &loadmx, &netmx, &frqmax, &angmax, 
               &mxcoup, &jmax, &nsmaxx, &npmax, &verbose, &debug, &remove_old );
  
  /* Verify input file name and file exists */
  job_name = (char*) malloc( strlen( nec_file ) + 2 );
  strncpy( job_name, nec_file, strlen( nec_file ) + 1 );
  verify_input_file( &job_name );

  /* Remove old output files */
  if( remove_old )
    remove_old_files( job_name );

  /* FNORM must be larger of angmax and frqmax */
  if( angmax >= frqmax )
    normf = angmax;
  else
    normf = frqmax;
  
  /* LD depends on number of segments and patches */
  d_num_seg = (double)num_seg;
  d_num_patch = (double)num_patch;
  d_ld = d_num_seg + d_num_patch;

  if( d_ld <= (double)MAX_LONG ) 
    {
      if ( (long)d_ld > LD_MAX ) 
	{
	  free( job_name );
	  error( 1, 0, "LD exceeds hard limit of %d\n", LD_MAX );
	}
      else 
	{
	  ld = (long)d_ld;
	}
    }
  else 
    {
      free( job_name );
      error( 1, 0, "LD too big for long\n" );
    }
  
  /* Verfiy core size and set IRESRV */
  in_core_iresrv = ( d_num_seg + 2 * d_num_patch ) * ( d_num_seg + 2 * d_num_patch );
  
  if( core_size == 0 ) 
    {
      d_iresrv = in_core_iresrv;
      core_size = 16.0 * d_iresrv / (double)ONE_MEGABYTE;
    }
  else 
    {
      d_iresrv = core_size * (double)ONE_MEGABYTE / 16.0;
      if( d_iresrv < 2.0 * ( (double)num_seg + 2.0 * (double)num_patch ) ) 
	{
	  free( job_name );
	  error( 1, 0, "core size too small for job" );
	}
      else if ( d_iresrv < 0.2 * in_core_iresrv ) 
	{
	  error( 0, 0, "warning: possibly massively out of core solution" );
	}
      else if ( d_iresrv < in_core_iresrv ) 
	{
	  error( 0, 0, "warning: possibly out of core solution" );
	}
    }
  
  if( d_iresrv <= (double)MAX_LONG ) 
    {
      iresrv = (long)d_iresrv;
    }
  else 
    {
      free( job_name );
      error( 1, 0, "core size too big for long" );
    }
  
  /* Dump values if verbose */
  if( verbose )
    dump_values( nec_file, multi_file, core_size, num_seg, num_patch, nsmax, 
		 loadmx, netmx, frqmax, angmax, mxcoup, jmax, nsmaxx, 
		 npmax, ld, iresrv, normf );
  
  /* Call NEC */
  ret = cnec2d( job_name, ld, iresrv, nsmax, loadmx, netmx, 
                normf, mxcoup, jmax , nsmaxx, npmax, multi_file, 
		verbose, debug );
  
  if( ret != CNEC2D_ERR_OK ) 
    {
      free( job_name );
      error( 1, 0, "%s", cnec2d_err( ret ) );
    }
  
  free( job_name );
  
  return 0;

}


/*
 * Give version info and exit.
 */

void 
version_info( int verbose )
{
  
  printf( "\ncnec2 (%s) %s\n", PACKAGE, VERSION );
  printf( "Copyright (C) 1999-2016 Ian David Flintoft\n" );
  printf( "Maintained by Ian D. Flintoft <ian.flintoft@googlemail.com>\n" );
  printf( "%s comes with NO WARRANTY,\n", PACKAGE );
  printf( "to the extent permitted by law.\n" );
  printf( "You may redistribute copies of %s\n", PACKAGE );
  printf( "under the terms of the GNU General Public License.\n" );
  printf( "For more information about these matters,\n" );
  printf( "see the file named COPYING.\n\n" );

#ifdef USE_BLAS
  printf( "Using LAPACK/BLAS LU solver.\n\n" );
#else
  printf( "Using NEC2 built in LU solver.\n\n" );
#endif
#ifdef BUILTIN_SOMNEC
  printf( "SOMNEC is built in.\n\n" );
#else
  printf( "SOMNEC support requires external program.\n\n" );
#endif

  printf( "Built on %s at %s.\n\n", __DATE__, __TIME__ );

  if( verbose )
    {
      printf( "-------------------------------------------------\n");
      printf( "      Item              Minimum  Default  Maximum\n");
      printf( "-------------------------------------------------\n");
      printf( "Number of segments     %8d %8d %8d\n", 
	      NUM_SEG_MIN, NUM_SEG_DEFAULT, NUM_SEG_MAX );
      printf( "Number of patches      %8d %8d %8d\n", 
	      NUM_PATCH_MIN, NUM_PATCH_DEFAULT, NUM_PATCH_MAX );
      printf( "Core size (Mb)         %8d %8d %8d\n", 
	      CORE_SIZE_MIN, CORE_SIZE_DEFAULT, CORE_SIZE_MAX );
      printf( "Number of sources      %8d %8d %8d\n", 
	      NSMAX_MIN, NSMAX_DEFAULT, NSMAX_MAX );
      printf( "Number of loads        %8d %8d %8d\n", 
	      LOADMX_MIN, LOADMX_DEFAULT, LOADMX_MAX );
      printf( "Number of networks     %8d %8d %8d\n", 
	      NETMX_MIN, NETMX_DEFAULT, NETMX_MAX );
      printf( "Number of frequencies  %8d %8d %8d\n", 
	      FRQMAX_MIN, FRQMAX_DEFAULT, FRQMAX_MAX );
      printf( "Number of angles       %8d %8d %8d\n", 
	      ANGMAX_MIN, ANGMAX_DEFAULT, ANGMAX_MAX );
      printf( "Number of couplings    %8d %8d %8d\n", 
	      MXCOUP_MIN, MXCOUP_DEFAULT, MXCOUP_MAX );
      printf( "Number of connections  %8d %8d %8d\n", 
	      JMAX_MIN, JMAX_DEFAULT, JMAX_MAX );
      printf( "Number of NGF segments %8d %8d %8d\n", 
	      NSMAXX_MIN, NSMAXX_DEFAULT, NSMAXX_MAX );
      printf( "Number of NGF patches  %8d %8d %8d\n", 
	      NPMAX_MIN, NPMAX_DEFAULT, NPMAX_MAX );
      printf( "-------------------------------------------------\n\n");
    }

  exit( EXIT_SUCCESS );

  return;

}


/* 
 * Give usage info and exit. 
 */

void 
usage( const int exitcode, const char *error, const char *add )
{
  
  printf( "\ncnec2 (%s) %s\n", PACKAGE, VERSION );
  printf( "\n" );
  printf( "Usage:  %s [options] jobname[.nec]\n\n", PACKAGE );
  printf( "where jobname.nec is the name of a NEC input file. The \".nec\" suffix\n" );
  printf( "may be omitted. The valid options are:\n\n" );
  printf( "    -h, --help            Give this help and exit\n" );
  printf( "    -V, --version         Give version info and exit\n" );
  printf( "    -v, --verbose         Give verbose output\n" );
  printf( "    -d, --debug=N         Debug level (default 0)\n" );
  printf( "    -s, --segments=N      Maximum number of segments is N (default 300)\n" );
  printf( "    -p, --patches=M       Maximum number of patches is M (default 0)\n" );
  printf( "    -c, --core-size=N     Limit core size to N bytes. The suffices b, k\n" );
  printf( "                          and m may be used for bytes, kilobytes and megabytes\n" );
  printf( "                          respectivley (default to exactly fit number of segments\n" );
  printf( "                          and patches specified with -s and -p options)\n" );
  printf( "    -e, --sources=N       Maximum number of sources is N (default 30)\n" );
  printf( "    -l, --loads=N         Maximum number of loads is N ( default 30)\n" );
  printf( "    -n, --networks=N      Maximum number of networks is N (default 10)\n" );
  printf( "    -x, --connect=N       Maximum number of segments meeting at a single\n" );  
  printf( "                          wire junction is N (default 30)\n" );
  printf( "    -u, --coupling=N      Maximum number of points in coupling calculation\n" );
  printf( "                          is N (default 5)\n" );
  printf( "    -S, --ngf-segments=N  Maximum number of NGF segments to which new\n" );
  printf( "                          segments or patches connect is N (default 50)\n" );
  printf( "    -P, --ngf-patches=N   Maximum number of NGF patches to which new\n" );
  printf( "                          segments or patches connect is N (default 50)\n" );
  printf( "    -f, --frequencies=N   Maximum number of frequencies in normalised\n" );
  printf( "                          impedance calculation is N (default 200)\n" );
  printf( "    -a, --angles=N        Maximum number angles in received signal\n" );
  printf( "                          strength is N (default 200)\n" );
  printf( "    -m, --multi-file      Enable multiple file output\n" );
  printf( "    -M, --old-multi-file  Enable old style multiple file output\n" );
  printf( "    -r, --remove-old      Remove old output files before running solver\n" );
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
get_options( int argc, char **argv , char **nec_file, int *multi_file, 
             double *core_size, long *num_seg, long *num_patch, 
             long *nsmax, long *loadmx, long *netmx, long *frqmax, 
	     long *angmax, long *mxcoup, long *jmax, long *nsmaxx, 
	     long *npmax, int *verbose, int *debug, int* remove_old )
{
  
  static struct option long_options[] = {
    { "version",           no_argument,       0, 'V' },
    { "help",              no_argument,       0, 'h' },
    { "verbose",           no_argument,       0, 'v' },
    { "debug",             required_argument, 0, 'd' },
    { "segments",          required_argument, 0, 's' },
    { "patches",           required_argument, 0, 'p' },
    { "core-size",         required_argument, 0, 'c' },
    { "sources",           required_argument, 0, 'e' },
    { "loads",             required_argument, 0, 'l' },
    { "networks",          required_argument, 0, 'n' },
    { "connect",           required_argument, 0, 'x' },
    { "coupling",          required_argument, 0, 'u' },
    { "ngf-segments",      required_argument, 0, 'S' },
    { "ngf-patches",       required_argument, 0, 'P' },
    { "frequencies",       required_argument, 0, 'f' },
    { "angles",            required_argument, 0, 'a' },
    { "multi-file",        no_argument,       0, 'm' },
    { "old-multi--file",   no_argument,       0, 'M' },
    { "remove-old",        no_argument,       0, 'r' },
    { 0, 0, 0, 0 }
  };
  char c;        /* Current short option */
  long tmp_long; /* Temporary long */

  while(1) 
    {
      
      c = getopt_long (argc, argv, "vhVd:s:p:c:e:l:n:x:u:S:P:f:a:mMr", long_options, 0 );
      
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
	case 's':
	  get_pos_long( c, num_seg, optarg, NUM_SEG_MIN, NUM_SEG_MAX );
	  break;
	case 'p':
	  get_pos_long( c, num_patch, optarg, NUM_PATCH_MIN, NUM_PATCH_MAX );
	  break;
	case 'c':
	  get_pos_bkm( c, core_size, optarg, CORE_SIZE_MIN, CORE_SIZE_MAX );
	  break;
	case 'e':
	  get_pos_long( c, nsmax, optarg, NSMAX_MIN, NSMAX_MAX );
	  break;
	case 'l':
	  get_pos_long( c, loadmx, optarg, LOADMX_MIN, LOADMX_MAX );
	  break;
	case 'n':
	  get_pos_long( c, netmx, optarg, NETMX_MIN, NETMX_MAX );
	  break;
	case 'x':
	  get_pos_long( c, jmax, optarg, JMAX_MIN, JMAX_MAX );
	  break;
	case 'u':
	  get_pos_long( c, mxcoup, optarg, MXCOUP_MIN, MXCOUP_MAX );
	  break;
	case 'S':
	  get_pos_long( c, nsmaxx, optarg, NSMAXX_MIN, NSMAXX_MAX );
	  break;
	case 'P':
	  get_pos_long( c, npmax, optarg, NPMAX_MIN, NPMAX_MAX );
	  break;
	case 'f':
	  get_pos_long( c, frqmax, optarg, FRQMAX_MIN, FRQMAX_MAX );
	  break;
	case 'a':
	  get_pos_long( c, angmax, optarg, ANGMAX_MIN, ANGMAX_MAX );
          break;
	case 'm':
	  *multi_file = 1;
	  break;
	case 'M':
	  *multi_file = 2;
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
  
  /* Get input file name as last arg */
  if( optind == argc )
    usage ( EXIT_FAILURE, "no input file name specified", "" );
  else if( optind == argc - 1 )
    *nec_file = argv[optind];
  else
    usage ( EXIT_FAILURE, "too many non-option arguments", "" );
  
  return;

}


/*
 * Parse a field consisting of a positive integer with an optional
 * suffix of b, k or m representing a sclae factor of 1, 1024 and
 * 1024*1024 respectively.
 *
 */

void
get_pos_bkm( const char c, double *val, const char *optarg, 
	     const long min, const long max )
{

  char *start;         /* Start of integer field */
  char *end;           /* End of integer field */
  char *junk;          /* Unused pointer in strtod */
  char *str;           /* Pointer to local copy of field */
  double scale = 1.0;  /* Scale factor from suffix */

  str = (char*) malloc( strlen( optarg ) + 1 );
  
  strncpy( str, optarg, strlen( optarg ) + 1 );

  start = str;

  if( *start == '+' ) start++;
  
  end = start;

  while( isdigit( *end ) ) end++;

  if( *end != '\0' ) 
    {
      switch( *end ) 
	{
	case 'b':
	  scale = 1.0;
	  *end = '\0';
      end++;
      break;
	case 'k':
	  scale = (double)ONE_KILOBYTE;
	  *end = '\0';
	  end++;
	  break;
	case 'm':
	  scale = (double)ONE_MEGABYTE;
	  *end = '\0';
	  end++;
	  break;
	default:
	  free( str );
	  error( 1, 0, "invalid value %s for option %c", optarg, c );
	  break;
	}
    }
  
  if( *end != '\0' ) 
    {
      free( str );
      error( 1, 0, "invalid value %s for option %c", optarg, c );
    }
  
  *val = strtod( start, &junk );
  *val *= scale;
  *val /= (double)ONE_MEGABYTE;
  
  if( *val < (long)min  || *val > (long)max ) 
    {
      free( str );
      error( 1, 0, "value of option %c must be in range %ldMb to %ldMb", 
	     c, min, max );
    }
  
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
 * Check input file exists and is openable.
 *
 * The name on the command line may have NEC_SUFFIX or not but
 * the actual file is required to have NEC_SUFFIX. Test for 
 * suffix on the file name - if not there add it.
 *
 */

void verify_input_file( char **job_name )
{

  FILE *fp;          /* Input file pointer */
  char *p;           /* Pointer to start of suffix */
  char *nec_file;    /* Pointer to file name */
  char *NEC_SUFFIX;  /* Input file suffix with dot */

  NEC_SUFFIX = (char *) malloc( strlen( SUFINP ) + 2 );
  strncpy( NEC_SUFFIX, ".", 2 );
  strncat( NEC_SUFFIX, SUFINP, strlen( SUFINP ) + 1 );
  p = strstr( *job_name, NEC_SUFFIX );

  if( p ) 
    {
      /* Suffix in file name */
      if( strncmp ( *job_name + strlen( *job_name ) - strlen( NEC_SUFFIX ), 
                    NEC_SUFFIX , strlen( NEC_SUFFIX ) + 1 ) != 0 ) 
	{
	  /* But not at end so ignore */
	  nec_file = (char *) malloc( strlen( *job_name ) + strlen( NEC_SUFFIX ) + 1 );
	  strncpy( nec_file, *job_name, strlen( *job_name ) + 1 );
	  strncat( nec_file, NEC_SUFFIX, strlen( NEC_SUFFIX ) + 1 );
	}
    else 
      {
	/* Strip off suffix for job name */
	nec_file = (char *) malloc( strlen( *job_name ) + 1 );
	strncpy( nec_file, *job_name, strlen( *job_name ) + 1 );
	*p = '\0';
      }
    }
  else 
    {
      /* No suffix - put one on to file name */
      nec_file = (char *) malloc( strlen( *job_name ) + strlen( NEC_SUFFIX ) + 1 );
      strncpy( nec_file, *job_name, strlen( *job_name ) + 1 );
      strncat( nec_file, NEC_SUFFIX, strlen( NEC_SUFFIX ) + 1 );
    }

  free( NEC_SUFFIX );

  /* Now try and open it */ 
  fp = fopen( nec_file, "r" );
  if( !fp ) 
    {
      free( nec_file );
      error( 1, 0, "cannot open input file %s", nec_file );    
    }

  if( ferror( fp ) )
    {
      free( nec_file );
      error( 1, 0, "error on input file %s", nec_file );    
    }

  if( fclose( fp ) == EOF )
    {
      free( nec_file );
      error( 1, 0, "cannot close input file %s", nec_file );    
    }
  
  return;
  
}


/* 
 * Dump out input values to stdout.
 *
 */

void 
dump_values( const char *nec_file, const int multi_file, const double core_size,
	     const long num_seg, const long num_patch, const long nsmax, 
	     const long loadmx, const long netmx, const long frqmax, 
	     const long angmax, const long mxcoup, const long jmax, 
	     const long nsmaxx, const long npmax, const long ld, 
	     const long iresrv, const long normf )
{

  printf( "\n\n" );
  printf( "Input file:             %s\n",  nec_file );
  printf( "Multi-file output:      %d\n",  multi_file );
  printf( "Core size (Mb):         %.1f\n", core_size );
  printf( "Number of segments:     %ld\n", num_seg );
  printf( "Number of patches:      %ld\n", num_patch );
  printf( "Number of sources:      %ld\n", nsmax );
  printf( "Number of loads:        %ld\n", loadmx );
  printf( "Number of networks:     %ld\n", netmx );
  printf( "Number of frequencies:  %ld\n", frqmax );
  printf( "Number of angles:       %ld\n", angmax );
  printf( "Value of NORMF:         %ld\n", normf );
  printf( "Number of couplings:    %ld\n", mxcoup );
  printf( "Number of connections:  %ld\n", jmax );
  printf( "Number of NGF segments: %ld\n", nsmaxx );
  printf( "Number of NGF patches:  %ld\n", npmax );
  printf( "Value of LD:            %ld\n", ld );
  printf( "Value of IRESRV:        %ld\n", iresrv );
  printf( "\n\n" );

  return;

}


/* 
 * Remove old output files.
 *
 */

void 
remove_old_files( const char *job_name )
{

  FILE *fpout;                     /* Output file pointer */
  char buf[10];                    /* Buffer for -m number */
  char *outfile  = 0;              /* Output file buffer */
  char suffices1[][5] = { ".res" , 
                          ".imp" , 
                          "" };    /* Suffices for -M files */
  char suffices2[][5] = { ".nre" , 
                          ".nrh" , 
                          ".rdp" , 
                          "" };    /* Suffice for -m files */
  int ret;                         /* Return code from remove */
  int i;                           /* Counter for suffix arrays */
  int j;                           /* Counter for -m files */
  int not_got_to_end;              /* True while more -m files left */
  
  outfile = (char *) malloc( strlen( job_name ) + CNEC2D_BUF_EXTRA );
  
  i=0;
  
  /* Delete files associated with -M option */
  
  while( suffices1[i][0] )
    {
      
      strncpy( outfile , job_name , strlen( job_name ) + 1 );
      strncat( outfile , suffices1[i] , strlen( suffices1[i] ) + 1 );  
      fpout = fopen( outfile , "r" );
      if ( fpout )
	{
	  fclose( fpout );
	  ret = remove( outfile );
	}
      
      i++;
    }
  
  i=0;
  
  while( suffices2[i][0] )
    {
      
      strncpy( outfile , job_name , strlen( job_name ) + 1 );
      strncat( outfile , suffices2[i] , strlen( suffices2[i] ) + 1 );  
      fpout = fopen( outfile , "r" );
      if ( fpout )
	{
	  fclose( fpout );
	  ret = remove( outfile );
	}
      
      i++;
    }
  
  i=0;
  
  /* Delete files associated with -m option */
  
  while( suffices2[i][0] )
    {
      
      j = 1;
      not_got_to_end = 1;
      
      do 
	{
	  
	  strncpy( outfile , job_name , strlen( job_name ) + 1 );
          sprintf( buf , "-%04d" , j );
	  strncat( outfile , buf , strlen( buf ) + 1 );  
	  strncat( outfile , suffices2[i] , strlen( suffices2[i] ) + 1 );  
	  fpout = fopen( outfile , "r" );
	  if ( fpout )
	    {
	      fclose( fpout );
	      ret = remove( outfile );
	    }
          else
	    not_got_to_end = 0;
	  
	  j++;
	  
	} while( not_got_to_end );
      
      i++;
    }
  
  return;
  
}
