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
 * cnec2d.c -- C driver routine for NEC2D.
 *
 * This file contains the functions for driving NEC2D with dynamically
 * allocated memory.
 *
 */

#include <config.h>
#include <fcint.h>
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cnec2d.h"

/* 
 * Prototype for FORTRAN routine assuming one underscore and length of character 
 * array is passed at end of arg list by value as an unsigned long. This is not
 * portable but works on many machines with many compilers.
 */
void FC_nec2d( char*, char*, INTEGER*, INTEGER*, INTEGER*,
              INTEGER*, INTEGER*, INTEGER*, INTEGER*, 
              INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, 
              INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, 
              INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, 
              INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, INTEGER*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, REAL_8*, 
              COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, 
              COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, INTEGER*, INTEGER*, 
	      INTEGER*, COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, COMPLEX_16*, 
	      COMPLEX_16*, FTNLEN, FTNLEN );


/* 
 * Private functions prototypes.
 */

static int alloc_integer( INTEGER**, const INTEGER ); 
static int alloc_real_8( REAL_8**, const INTEGER ); 
static int alloc_complex_16( COMPLEX_16**, const INTEGER ); 


/*
 * Driver routine.
 */

int 
cnec2d( char *JOBNAM, const long ld, const long iresrv, const long nsmax, 
	const long loadmx, const long netmx, const long normf, const long mxcoup, 
	const long jmax , const long nsmaxx, const long npmax, const int multi_file,
	int verbose, int debug )
{

  char *FILBUF;
  FTNLEN jobnam_len, filbuf_len;
  INTEGER LD , IRESRV , NSMAX , LOADMX , NETMX , NORMF , 
          MXCOUP , MXCOU2 , JMAX , NSMAXX , NPMAX, MULTIF, 
          IFAIL , DEBUG;
  INTEGER *IP , *ICON1 , *ICON2 , *ITAG , *ICONX , *IVQD ,
          *ISANT , *IQDS , *NCTAG , *NCSEG , *ISEG1 , *ISEG2 ,
          *NTYP , *JCO , *ISCON , *IPCON , *LDTYP , *LDTAG , 
          *LDTAGF , *LDTAGT , *IX , *IPNT , *NTEQA , *NTSCA;
  REAL_8 *SALP , *AIR , *AII , *BIR , *BII , *CIR , *CII , *GAIN , 
         *X , *Y , *Z , *SI , *BI , *ALP , *BET , *X2 , *Y2 , *Z2 , 
         *CAB , *SAB , *T1X , *T1Y , *T1Z , *T2X , *T2Y , *T2Z , *XS , 
         *YS , *ZS , *Z2S , *X11R , *X11I , *X12R , *X12I , *X22R , 
         *X22I , *AX , *BX , *CX , *ZLR , *ZLI , *ZLC , *FNORM , *XTEMP , 
         *YTEMP , *ZTEMP , *SITEMP , *BITEMP;
  COMPLEX_16 *CM , *CUR , *ZARRAY , *D , *VQD , *VSANT , *VQDS , 
             *Y11A , *Y12A , *CMN , *RHNT , *RHS , *VSRC , *RHNX; 

  LD = (INTEGER) ld;
  IRESRV = (INTEGER) iresrv;
  NSMAX  = (INTEGER) nsmax;
  LOADMX = (INTEGER) loadmx;
  NETMX  = (INTEGER) netmx;
  NORMF  = (INTEGER) normf;
  MXCOUP = (INTEGER) mxcoup;
  MXCOU2 = MXCOUP * MXCOUP - MXCOUP;
  JMAX   = (INTEGER) jmax;
  NSMAXX = (INTEGER) nsmaxx;
  NPMAX  = (INTEGER) npmax;
  MULTIF = (INTEGER) multi_file;
  DEBUG  = (INTEGER) debug;

  /* Length of job name needed for call to FORTRAN routine. */
  jobnam_len = (FTNLEN) strlen( JOBNAM );

  /* Filename buffer for FORTRAN routines. CNEC2D_BUF_EXTRA is */
  /* at least size of suffix + number digits used by FILSUF plus 3. */
  /* CAREFUL. Fortran must be passed the size of the buffer NOT the */
  /* length of the string in it, which is zero! */
  FILBUF = (char *) calloc( jobnam_len + (FTNLEN)CNEC2D_BUF_EXTRA , sizeof(char) );
  filbuf_len = jobnam_len + (FTNLEN)CNEC2D_BUF_EXTRA;

  /* Allocate NEC arrays */
  alloc_integer( &IP, 2 * LD );
  alloc_integer( &ICONX, 2 * LD );
  alloc_integer( &IVQD, NSMAX );
  alloc_integer( &ISANT, NSMAX );
  alloc_integer( &IQDS, NSMAX );
  alloc_integer( &NCTAG, MXCOUP );
  alloc_integer( &NCSEG, MXCOUP );
  alloc_integer( &ISEG1, NETMX );
  alloc_integer( &ISEG2, NETMX );
  alloc_integer( &NTYP, NETMX );
  alloc_integer( &JCO, JMAX );
  alloc_integer( &ISCON, NSMAXX );
  alloc_integer( &IPCON, NPMAX );
  alloc_integer( &LDTYP, LOADMX );
  alloc_integer( &LDTAG, LOADMX );
  alloc_integer( &LDTAGF, LOADMX );
  alloc_integer( &LDTAGT, LOADMX );
  alloc_integer( &IX, 2 * LD );
  alloc_real_8( &SALP, LD );
  alloc_real_8( &AIR, LD );
  alloc_real_8( &AII, LD );
  alloc_real_8( &BIR, LD );
  alloc_real_8( &BII, LD );
  alloc_real_8( &CIR, LD );
  alloc_real_8( &CII, LD );
  alloc_real_8( &GAIN, 4 * LD );
  alloc_real_8( &X11R, NETMX );
  alloc_real_8( &X11I, NETMX );
  alloc_real_8( &X12R, NETMX );
  alloc_real_8( &X12I, NETMX );
  alloc_real_8( &X22R, NETMX );
  alloc_real_8( &X22I, NETMX );
  alloc_real_8( &AX, JMAX );
  alloc_real_8( &BX, JMAX );
  alloc_real_8( &CX, JMAX );
  alloc_real_8( &ZLR, LOADMX );
  alloc_real_8( &ZLI, LOADMX );
  alloc_real_8( &ZLC, LOADMX );
  alloc_real_8( &FNORM, NORMF );
  alloc_real_8( &XTEMP, LD );
  alloc_real_8( &YTEMP, LD );
  alloc_real_8( &ZTEMP, LD );
  alloc_real_8( &SITEMP, LD );
  alloc_real_8( &BITEMP, LD );
  alloc_complex_16( &CM, IRESRV );
  alloc_complex_16( &CUR, 3 * LD );
  alloc_complex_16( &ZARRAY, LD );
  alloc_complex_16( &D, 2 * LD );
  alloc_complex_16( &VQD, NSMAX );
  alloc_complex_16( &VSANT, NSMAX );
  alloc_complex_16( &VQDS, NSMAX );
  alloc_complex_16( &Y11A, MXCOUP );
  alloc_complex_16( &Y12A, MXCOU2 );
  alloc_integer( &IPNT, NETMX );
  alloc_integer( &NTEQA, NETMX );
  alloc_integer( &NTSCA, NETMX );
  alloc_complex_16( &CMN, NETMX * NETMX );
  alloc_complex_16( &RHNT, NETMX );
  alloc_complex_16( &RHS, 3 * LD );
  alloc_complex_16( &VSRC, NETMX );
  alloc_complex_16( &RHNX, NETMX );

  /* Dynamic arrays with equivalent storage */
  alloc_real_8( &SI, LD );
  X2 = SI;
  T1X = SI;

  alloc_real_8( &ALP, LD );
  Y2 = ALP;
  CAB = ALP;
  T1Y = ALP;

  alloc_real_8( &BET, LD );
  Z2 = BET;
  SAB = BET;
  T1Z = BET;

  alloc_integer( &ICON1, 2 * LD );
  T2X = (REAL_8 *) ICON1;

  alloc_integer( &ICON2, 2 * LD );
  T2Y = (REAL_8 *) ICON2;

  alloc_integer( &ITAG, 2 * LD );
  T2Z = (REAL_8 *) ITAG;

  alloc_real_8( &X, LD );
  XS = X;

  alloc_real_8( &Y, LD );
  YS = Y;

  alloc_real_8( &Z, LD );
  ZS = Z;

  alloc_real_8( &BI, LD );
  Z2S = BI;

  /* Call NEC2 */
  IFAIL = 0;
  FC_nec2d(JOBNAM,FILBUF,&IFAIL,&DEBUG,&MULTIF,&LD,&IRESRV,&NSMAX,&LOADMX,
            &NETMX,&NORMF,&MXCOUP, &MXCOU2,&JMAX,&NSMAXX,&NPMAX,
           IP,ICON1,ICON2,ITAG,ICONX,IVQD,ISANT,IQDS,NCTAG,NCSEG,
           ISEG1,ISEG2,NTYP,JCO,ISCON,IPCON,LDTYP,LDTAG,LDTAGF,
           LDTAGT,IX,SALP,AIR,AII,BIR,BII,CIR,CII,GAIN,X,Y,Z,SI,
           BI,ALP,BET,X2,Y2,Z2,CAB,SAB,T1X,T1Y,T1Z,T2X,T2Y,T2Z,
           XS,YS,ZS,Z2S,X11R,X11I,X12R,X12I,X22R,X22I,AX,BX,CX,
           ZLR,ZLI,ZLC,FNORM,XTEMP,YTEMP,ZTEMP,SITEMP,BITEMP,CM,
           CUR,ZARRAY,D,VQD,VSANT,VQDS,Y11A,Y12A,IPNT,NTEQA,NTSCA,
	   CMN,RHNT,RHS,VSRC,RHNX,jobnam_len,filbuf_len);

  /* Free up memory */
  free( IP );
  free( ICONX );
  free( IVQD );
  free( ISANT );
  free( IQDS );
  free( NCTAG );
  free( NCSEG );
  free( ISEG1 );
  free( ISEG2 );
  free( NTYP );
  free( JCO );
  free( ISCON );
  free( IPCON );
  free( LDTYP );
  free( LDTAG );
  free( LDTAGF );
  free( LDTAGT );
  free( IX );
  free( SALP );
  free( AIR );
  free( AII );
  free( BIR );
  free( BII );
  free( CIR );
  free( CII );
  free( GAIN );
  free( X11R );
  free( X11I );
  free( X12R );
  free( X12I );
  free( X22R );
  free( X22I );
  free( AX );
  free( BX );
  free( CX );
  free( ZLR );
  free( ZLI );
  free( ZLC );
  free( FNORM );
  free( XTEMP );
  free( YTEMP );
  free( ZTEMP );
  free( SITEMP );
  free( BITEMP );
  free( CM );
  free( CUR );
  free( ZARRAY );
  free( D );
  free( VQD );
  free( VSANT );
  free( VQDS );
  free( Y11A );
  free( Y12A );
  free( IPNT );
  free( NTEQA );
  free( NTSCA );
  free( CMN );
  free( RHNT );
  free( RHS );
  free( VSRC );
  free( RHNX );
  free( SI );
  free( ALP );
  free( BET );
  free( ICON1 );
  free( ICON2 );
  free( ITAG );
  free( X );
  free( Y );
  free( Z );
  free( BI );
  free( FILBUF );

  return IFAIL;

}


/*
 * Allocate an INTEGER array.
 */

static int 
alloc_integer( INTEGER **a, const INTEGER size )
{

  *a = (INTEGER *) calloc( size , sizeof(INTEGER) );

  return 0;

}


/*
 * Allocate a REAL_8 array.
 */

static int 
alloc_real_8( REAL_8 **a, const INTEGER size )
{

  *a = (REAL_8 *) calloc( size , sizeof(REAL_8) );

  return 0;

}


/*
 * Allocate a COMPLEX_16 array.
 */

static int 
alloc_complex_16( COMPLEX_16 **a, const INTEGER size )
{
  
  *a = (COMPLEX_16 *) calloc( size , sizeof(COMPLEX_16) );
  
  return 0;
  
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

static char err_msg[CNEC2D_ERR_END][MAX_ERR_LEN] = {
  "successful execution",                                                                                                                                                /* IFAIL =  0 */
  "malloc error in cnec2d",                                                                                                                                              /* IFAIL =  1 */
  "arc angle exceeds 360 degrees",                                                                                                                                       /* IFAIL =  2 */
  "error reading matrix from file",                                                                                                                                      /* IFAIL =  3 */
  "range error in conect\n\nAn end of file has been encountered while reading data from the\nmatrix file.",                                                             /* IFAIL =  4 */
  "segment extends below ground\n\nWhen ground is specified on the GE card, no segment may extend\nbelow the XY plane.",                                                /* IFAIL =  5 */
  "segment lies in ground plane\n\nWhen ground is specified on the GE card, no segment should lie\nin the XY plane.",                                                   /* IFAIL =  6 */
  "number of new segments or patches in NGF solution exceeds limit",                                                                                                     /* IFAIL =  7 */
  "invalid value for ipsym in conect",                                                                                                                                   /* IFAIL =  8 */
  "segment connection error\n\nPossible causes: number of segments at a junction exceeds limit;\nsegment lengths are zero; array overflow.",			         /* IFAIL =  9 */
  "geometry data card error\n\nA geometry data card was expected, but the card mnemonic is not\nthat of a geometry card; check GE present.",			         /* IFAIL = 10 */
  "patch data error\n\nInvalid data on SP, SM, or SC card; or SC card not found where required.",							                 /* IFAIL = 11 */
  "GF must be first geometry data card",				                                                                                                 /* IFAIL = 12 */
  "segment data error\n\nA segment with zero length or zero radius was found.",				                                                                 /* IFAIL = 13 */
  "number of unknowns exceeds limit\n\nThe sum of the number of segments and patches exceeds the allowed\nnumber of unknowns.",	                                 /* IFAIL = 14 */
  "error in fblock",                                                                                                                                                     /* IFAIL = 15 */
  "insufficient storage for matrix\n\nArray storage for matrix is not sufficient for out-of-core solution.",                                                             /* IFAIL = 16 */
  "symmetry error\n\nArray overflow or program malfunction.",							                                                         /* IFAIL = 17 */
  "insufficient storage for interaction matrices\n\nArray storage exceeded in NGF solution.",			                                                         /* IFAIL = 18 */
  "internal error: file name too long for buffer",			                                                                                                 /* IFAIL = 19 */
  "too many output files",						                                                                                                 /* IFAIL = 20 */
  "internal error: file number less than zero",				                                                                                                 /* IFAIL = 21 */
  "parameter specifying segment position in a group of equal tags cannot be zero\n\nMay occur at any point where a tag number is used to identify a segment.",           /* IFAIL = 22 */
  "invalid segment number\n\nThis error results from faulty input data and can occur at any point\nwhere a tag number is used to identify a segment.",                  /* IFAIL = 23 */
  "improper load type\n\nValid load types (LDTYP on the LD card) are from 0 through 5.",                                                                                 /* IFAIL = 24 */
  "loading may not be added to segments in NGF section",		                                                                                                 /* IFAIL = 25 */
  "invalid segment number in loading card\n\nITAG specified on an LD card could not be found as a segment tag.",                                                         /* IFAIL = 26 */
  "network array dimensions too small\n\nThe number of different segments to which transmission lines or network\nports are connected exceeds array dimensions.",       /* IFAIL = 27 */
  "parse error",							                                                                                                 /* IFAIL = 28 */
  "four corners of quadrilateral patch do not lie in a plane\n\nThe four corners of a quadrilateral patch (SP card) must lie in a plane.",		                 /* IFAIL = 29 */
  "segment lies in plane of symmetry\n\nA segment may not lie in or cross a plane of symmetry about which the\nstructure is reflected.",                                /* IFAIL = 30 */
  "patch lies in plane of symmetry",		                                                                                            			         /* IFAIL = 31 */
  "B less than A in rom2",                                                                                                                                               /* IFAIL = 32 */
  "segment connection error\n\nThe number of segments at a junction exceeds dimension limit, or the\nconnection numbers are not self-consistant.",                      /* IFAIL = 33 */
  "error opening input file",						                                                                                                 /* IFAIL = 34 */
  "incorrect label for comment card\n\nComment cards must be the first cards in a data set, and the comments must\nbe terminated by the CE mnemonic.",		         /* IFAIL = 35 */
  "invalid data card after geometry section\n\nA card with an unrecognizable mnemonic was found in the program\ncontrol cards following the geometry cards.",		 /* IFAIL = 36 */
  "FR/GN card not allowed for NGF solution",				                                                                                                 /* IFAIL = 37 */
  "number of segments in coupling calculation exceeds limit",		                                                                                                 /* IFAIL = 38 */
  "number of loading cards exceeds limit\n\nThe number of LD cards exceeds array dimension.",                                                                            /* IFAIL = 39 */
  "data fault in loading card\n\nWhen several segments are loaded, the number of the second segment\nspecified must be greater than the number of the first segment.",  /* IFAIL = 40 */
  "radial wire GS approx may not be used with Sommerfeld ground", 	                                                                                                 /* IFAIL = 41 */
  "number of excitation cards exceeds limit\n\nThe number of voltage source excitations exceeds array dimensions.",	                                                 /* IFAIL = 42 */
  "number of network cards exceeds limit\n\nThe number of NT and TL cards exceeds array dimension.",	                                                                 /* IFAIL = 43 */
  "NGF already in use, cannot write new NGF",			                                                         		                                 /* IFAIL = 44 */
  "error opening Sommerfeld ground file",				                                                                                                 /* IFAIL = 45 */
  "error in ground parameters\n\nComplex dielectric constant from somnec file does not agree with\ndata from GN and FR cards.",	                                 /* IFAIL = 46 */
  "too many unknowns generated while executing GM card",			                                                                                         /* IFAIL = 47 */
  "NGF file already exists - please delete and rerun",			                                                                                                 /* IFAIL = 48 */
  "Results file already exists - please delete and rerun"			                                                                                         /* IFAIL = 49 */
};

/*
 * Return pointer to string containing error message for a given
 * error code.
 */

char *
cnec2d_err( const int err_code )
{
  
  if( err_code < CNEC2D_ERR_END )
    {
      return err_msg[err_code];
    }
  else 
    {
      fprintf( stderr, "*** Internal error: error code too big - please report bug! ***" );
      exit( EXIT_FAILURE );
    }
  
}
