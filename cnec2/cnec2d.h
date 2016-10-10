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
 * cnec2d.h -- C driver interface for NEC2D.
 *
 * This file contains the interface definition for the NEC2D C driver.
 *
 */

/* Mapping C and Fortran types */

typedef int INTEGER;
typedef double REAL_8;
struct COMPLEX_16_TAG {
  REAL_8 re;
  REAL_8 im;
};
typedef struct COMPLEX_16_TAG COMPLEX_16;
typedef unsigned long FTNLEN;

/* 
 * Return codes. Corresponding messages are in the static array in
 * cnec2d.c in the same order.
 */

enum {
  CNEC2D_ERR_OK = 0,                 /*  0 Successful execution */
  CNEC2D_ERR_MALLOC,                 /*  1 Malloc error in cnec2d */
  CNEC2D_ERR_ARCANG,                 /*  2 Arc angle exceeds 360 degrees */
  CNEC2D_ERR_BLCKIN,                 /*  3 Error reading matrix from file */
  CNEC2D_ERR_CONECT,                 /*  4 Range error in conect */
  CNEC2D_ERR_SEGGND,                 /*  5 Segment extends below ground */
  CNEC2D_ERR_SEGGPL,                 /*  6 Segment lies in ground plane */
  CNEC2D_ERR_SEGNGF,                 /*  7 Number of new segments connected to NGF segments or patches exceeds limit */
  CNEC2D_ERR_CONSYM,                 /*  8 Invalid value for ipsym in conect */
  CNEC2D_ERR_CONERR,                 /*  9 Segment connection error */
  CNEC2D_ERR_GEOERR,                 /* 10 Geometry data card error */
  CNEC2D_ERR_PATERR,                 /* 11 Patch data error */
  CNEC2D_ERR_GFORDR,                 /* 12 GF must be first geometry data card */
  CNEC2D_ERR_SEGERR,                 /* 13 Segment data error */
  CNEC2D_ERR_SEGLIM,                 /* 14 Number of unknowns exceeds limit */
  CNEC2D_ERR_FBLOCK,                 /* 15 Error in fblock */
  CNEC2D_ERR_CORELM,                 /* 16 Insufficient storage for matrix */
  CNEC2D_ERR_SYMERR,                 /* 17 Symmetry error */
  CNEC2D_ERR_INTLIM,                 /* 18 Insufficient storage for interaction matrices */
  CNEC2D_ERR_FILEMP,                 /* 19 Empty file name */
  CNEC2D_ERR_FILSUF,                 /* 20 No suffix on file name */
  CNEC2D_ERR_FILLNG,                 /* 21 New file name suffix too long */
  CNEC2D_ERR_SEGPAR,                 /* 22 Parameter specifying segment position in a group of equal tags must not be zero */
  CNEC2D_ERR_INVSEG,                 /* 23 Invalid segment number */
  CNEC2D_ERR_LODTYP,                 /* 24 Improper load type */
  CNEC2D_ERR_LODNGF,                 /* 25 Loading may not be added to segments in NGF section */
  CNEC2D_ERR_LODERR,                 /* 26 Invalid segment number in loading card */
  CNEC2D_ERR_NETLIM,                 /* 27 Network array dimensions too small */
  CNEC2D_ERR_CARDER,                 /* 28 Parse error */
  CNEC2D_ERR_PATQUA,                 /* 29 Corners of quadrilateral patch do not lie in a plane */
  CNEC2D_ERR_SYMSEG,                 /* 30 Segment lies in plane of symmetry */
  CNEC2D_ERR_SYMPAT,                 /* 31 Patch lies in plane of symmetry */
  CNEC2D_ERR_ROMBEA,                 /* 32 B less than A in rom2 */
  CNEC2D_ERR_SEGCON,                 /* 33 Segment connection error */
  CNEC2D_ERR_FILINP,                 /* 34 Error opening input file */
  CNEC2D_ERR_COMLAB,                 /* 35 Incorrect label for comment card */
  CNEC2D_ERR_GEOLAB,                 /* 36 Fautly data card label after geometry section */
  CNEC2D_ERR_CRDNGF,                 /* 37 Invalid card for NGF solution */
  CNEC2D_ERR_SEGCOU,                 /* 38 Number of segments in coupling calculation exceeds limit */
  CNEC2D_ERR_LODLIM,                 /* 39 Number of loading cards exceeds limit */
  CNEC2D_ERR_LODDAT,                 /* 40 Data fault in loading card */
  CNEC2D_ERR_RADSOM,                 /* 41 Radial wire GS approx may not be used with Sommerfeld ground */ 
  CNEC2D_ERR_EXTLIM,                 /* 42 Number of excitation cards exceeds limit */
  CNEC2D_ERR_NTCLIM,                 /* 43 Number of network cards exceeds limit */
  CNEC2D_ERR_NGFUSE,                 /* 44 NGF already in use */
  CNEC2D_ERR_SOMOPN,                 /* 45 Error opening Sommerfeld ground file */
  CNEC2D_ERR_SOMGND,                 /* 46 Error in ground parameters */
  CNEC2D_ERR_MOVLIM,                 /* 47 Too many unknowns generated while executing GM card */
  CNEC2D_ERR_NGFEXT,                 /* 48 NGF file already exists - please delete and rerun */
  CNEC2D_ERR_RESEXT,                 /* 49 Results file already exists - please delete and rerun */
  CNEC2D_ERR_END                     /* End marker */
};

/*
 * Extra length of file names used by NEC2D at the end of JOBNAM. Should
 * be at least size of suffix + number digits used by FILSUF + 3 ('.', '-' 
 * and '\0'). Currently this is 3 + 4 + 3 = 10. 
 */

enum {
  CNEC2D_BUF_EXTRA = 15
};

/* 
 * External interface.
 */

int cnec2d( char *, const long, const long, const long , const long, 
            const long, const long, const long, const long, const long, 
	    const long, const int, const int, const int );

char *cnec2d_err( const int );

