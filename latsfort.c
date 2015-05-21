/* -*-Mode: C;-*-
 * Module:      LATS Fortran interface
 *
 * Copyright:	1996, Regents of the University of California
 *		This software may not be distributed to others without
 *		permission of the author.
 *
 * Author:      Bob Drach, Lawrence Livermore National Laboratory
 *              drach@llnl.gov
 *
 * Version:     $Id: latsfort.c,v 1.10 1996/12/12 18:39:44 fiorino Exp $
 *
 * Revision History:
 *
 * $Log: latsfort.c,v $
 * Revision 1.10  1996/12/12 18:39:44  fiorino
 * 961212
 *
 * Mike's changes to parm table
 * GraDS updated source
 * and improvements in GrADS output to handle yearly data and 365-day calendars
 * added gaprnt routine in latsfort.c to aid link with straight GraDS source
 *
 * Revision 1.9  1996/11/11 22:39:19  drach
 * - Added function to set the basetime (lats_basetime)
 *
 * Revision 1.8  1996/10/22 19:05:06  fiorino
 * latsgrib bug in .ctl creator
 *
 * Revision 1.7  1996/08/27 19:41:39  drach
 * - Handle Cray double/double precision mismatch
 *
 * Revision 1.6  1996/08/20 18:34:08  drach
 * - lats_create has a new argument: calendar
 * - lats_grid: longitude, latitude dimension vectors are now double
 *   precision (double, C).
 * - lats_vert_dim: vector of levels is now double precision (double,
 *   C). lats_vert_dim need not be called for single-value/surface
 *   dimensions, if defined in the parameter table. Multi-valued vertical
 *   dimensions, such as pressure levels, must be defined with
 *   lats_vert_dim.
 * - lats_var: set level ID to 0 for implicitly defined surface
 *   dimension.
 * - lats_write: level value is double precision (double, C).
 * - lats_parmtab: replaces routine lats_vartab.
 * - FORTRAN latserropt added: allows program to set error handling
 *   options.
 * - The parameter file format changed.
 *
 * Revision 1.5  1996/06/27 01:11:10  drach
 * - Check for POSIX compliance
 *
 * Revision 1.4  1996/05/25 00:27:49  drach
 * - Added tables for vertical dimension types, time statistics, originating
 *   centers, and quality control marks
 * - Modified signatures of lats_create and lats_vert_dim
 *
 * Revision 1.3  1996/05/10 22:44:40  drach
 * - Initial version before GRIB driver added:
 * - Made grids, vertical dimensions file-independent
 *
 * Revision 1.2  1996/05/04 01:11:11  drach
 * - Added name, units to lats_vert_dim
 * - Added missing data attribute (latsnc.c)
 *
 * Revision 1.1  1996/05/03 18:48:06  drach
 * - Initial repository version
 *
 *
 */

#define _POSIX_SOURCE 1
#include "latsint.h"
#include "cfortran.h"

FCALLSCFUN5(INT,lats_basetime,LATSBASETIME,latsbasetime,INT,INT,INT,INT,INT)

FCALLSCFUN1(INT,lats_close,LATSCLOSE,latsclose,INT)

FCALLSCFUN8(INT,lats_create,LATSCREATE,latscreate,STRING,INT,INT,INT,INT,STRING,STRING,STRING)

FCALLSCSUB1(lats_erropt, LATSERROPT, latserropt, INT)

#ifdef CRAY
					     /* Handle double/DOUBLE PRECISION Cray mismatch */
int lats_grid_ld(char *name, int gridtype, int nlon, long double lons[], int nlat, long double lats[]){
	extern char *_lats_routine_name_;
	double *dlons, *dlats, *d;
	int i;

	_lats_routine_name_ = "lats_grid";

	if((dlons = (double *)calloc(nlon, sizeof(double)))==0 ||
	   (dlats = (double *)calloc(nlat, sizeof(double)))==0){
		latsError("Allocating lon/lat vectors, size %d, %d", nlon, nlat);
		return 0;
	}

					     /* Copy long double vectors to doubles */
	for(i=0, d=dlons; i<nlon; i++){
		*d++ = (double)lons[i];
	}

	for(i=0, d=dlats; i<nlat; i++){
		*d++ = (double)lats[i];
	}
	if(lats_grid(name, gridtype, nlon, dlons, nlat, dlats)==0)
		return 0;

	free(dlons);
	free(dlats);
	return 1;
}


FCALLSCFUN6(INT,lats_grid_ld,LATSGRID,latsgrid,STRING,INT,INT,DOUBLEV,INT,DOUBLEV)
#else
FCALLSCFUN6(INT,lats_grid,LATSGRID,latsgrid,STRING,INT,INT,DOUBLEV,INT,DOUBLEV)
#endif

FCALLSCFUN4(INT,lats_miss_float,LATSMISSREAL,latsmissreal,INT,INT,FLOAT,FLOAT)

FCALLSCFUN3(INT,lats_miss_int,LATSMISSINT,latsmissint,INT,INT,INT)

FCALLSCFUN7(INT,lats_var,LATSVAR,latsvar,INT,STRING,INT,INT,INT,INT,STRING)

FCALLSCFUN1(INT,lats_parmtab,LATSPARMTAB,latsparmtab,STRING)

FCALLSCSUB1(lats_qcopt,LATSQCOPT,latsqcopt,INT)

#ifdef CRAY
					     /* Handle double/DOUBLE PRECISION Cray mismatch */
int lats_vert_dim_ld(char* name, char* type, int nlev, long double levs[]){
	extern char *_lats_routine_name_;

	double *dlevs, *d;
	int i;

	if((dlevs = (double *)calloc(nlev, sizeof(double)))==0){
		latsError("Allocating lev vectors, size %d", nlev);
		return 0;
	}
					     /* Copy long double vectors to doubles */
	for(i=0, d=dlevs; i<nlev; i++){
		*d++ = (double)levs[i];
	}

	if(lats_vert_dim(name, type, nlev, dlevs)==0)
		return 0;

	free(dlevs);
	return 1;
}

FCALLSCFUN4(INT,lats_vert_dim_ld,LATSVERTDIM,latsvertdim,STRING,STRING,INT,DOUBLEV)
#else
FCALLSCFUN4(INT,lats_vert_dim,LATSVERTDIM,latsvertdim,STRING,STRING,INT,DOUBLEV)
#endif

#ifdef CRAY
					     /* Handle double/DOUBLE PRECISION Cray mismatch */
int lats_write_ld(int fileid, int varid, long double lev, int year, int month, int day, int hour, void* data){
	return lats_write(fileid, varid, (double)lev, year, month, day, hour, data);
}
FCALLSCFUN8(INT,lats_write_ld,LATSWRITE,latswrite,INT,INT,DOUBLE,INT,INT,INT,INT,PVOID)
#else
FCALLSCFUN8(INT,lats_write,LATSWRITE,latswrite,INT,INT,DOUBLE,INT,INT,INT,INT,PVOID)
#endif

/*---
  dump this routine defined in gauser.c but needed for linking by the GrADS routines
---*/

/* jfp gaprnt, gxgsym deleted, also defined in libcdms gaprnt.c

void gaprnt (int i, char *ch) {
  printf ("%s",ch);
}
*/
/* Query env symbol */
/* jfp gaprnt, gxgsym deleted, also defined in libcdms gaprnt.c
char *gxgsym(char *ch) {
  return (getenv(ch));
}
*/

