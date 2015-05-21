#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
/* version 1.2.1 of grib headers  w. ebisuzaki */

#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 0x7f) << 8) + b))
#endif

#define BDS_LEN(bds)		((int) ((bds[0]<<16)+(bds[1]<<8)+bds[2]))
#define BDS_Flag(bds)		(bds[3])

#define BDS_Grid(bds)		((bds[3] & 128) == 0)
#define BDS_Harmonic(bds)	(bds[3] & 128)

#define BDS_Packing(bds)	((bds[3] & 64) != 0)
#define BDS_SimplePacking(bds)	((bds[3] & 64) == 0)
#define BDS_ComplexPacking(bds)	((bds[3] & 64) != 0)

#define BDS_OriginalType(bds)	((bds[3] & 32) != 0)
#define BDS_OriginalFloat(bds)	((bds[3] & 32) == 0)
#define BDS_OriginalInt(bds)	((bds[3] & 32) != 0)

#define BDS_MoreFlags(bds)      ((bds[3] & 16) != 0)
#define BDS_UnusedBits(bds)	((int) (bds[3] & 15))

#define BDS_BinScale(bds)	INT2(bds[4],bds[5])

#define BDS_RefValue(bds)	(ibm2flt(bds+6))
#define BDS_NumBits(bds)	((int) bds[10])

#define BDS_DataStart(bds)      ((int) (11 + BDS_MoreFlags(bds)*3))

/* breaks if BDS_NumBits(bds) == 0 */
   #define BDS_NValues(bds)        (((BDS_LEN(bds) - BDS_DataStart(bds))*8 - \
				BDS_UnusedBits(bds)) / BDS_NumBits(bds))
/*
#define BDS_NValues(bds)        ((BDS_NumBits(bds) == 0) ? 0 : \
				(((BDS_LEN(bds) - BDS_DataStart(bds))*8 - \
				BDS_UnusedBits(bds)) / BDS_NumBits(bds)))
*/

/* undefined value -- if bitmap */
#define UNDEFINED		9.999e20

/* version 1.2 of grib headers  w. ebisuzaki */

#define BMS_LEN(bms)		((bms) == NULL ? 0 : (bms[0]<<16)+(bms[1]<<8)+bms[2])
#define BMS_UnusedBits(bms)	((bms) == NULL ? 0 : bms[3])
#define BMS_StdMap(bms)		((bms) == NULL ? 0 : ((bms[4]<<8) + bms[5]))
#define BMS_bitmap(bms)		((bms) == NULL ? NULL : (bms)+6)
#define BMS_nxny(bms)		((((bms) == NULL) || BMS_StdMap(bms)) \
	? 0 : (BMS_LEN(bms)*8 - 48 - BMS_UnusedBits(bms)))
/* cnames_file.c */

/* cnames.c */
/* then default values */
char *k5toa_def(int i, int center);
int atok5_def(char *name, int center);

/* version 1.3 of grib headers  w. ebisuzaki */
/* this version is incomplete */
#ifndef INT3
#define INT3(a,b,c) ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 127) << 16)+(b<<8)+c))
#endif
#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 127) << 8) + b))
#endif

#define GDS_Len1(gds)		(gds[0])
#define GDS_Len2(gds)		(gds[1])
#define GDS_Len3(gds)		(gds[2])
#define GDS_LEN(gds)		((int) ((gds[0]<<16)+(gds[1]<<8)+gds[2]))

#define GDS_NV(gds)		(gds[3])
#define GDS_PV(gds)		(gds[4])
#define GDS_PL(gds)		(gds[4])
#define GDS_DataType(gds)	(gds[5])

#define GDS_LatLon(gds)		(gds[5] == 0)
#define GDS_Mercator(gds)	(gds[5] == 1)
#define GDS_Gnomonic(gds)	(gds[5] == 2)
#define GDS_Lambert(gds)	(gds[5] == 3)
#define GDS_Gaussian(gds)	(gds[5] == 4)
#define GDS_Polar(gds)		(gds[5] == 5)
#define GDS_Harmonic(gds)	(gds[5] == 50)
#define GDS_Generic(gds)	(gds[5] == 220) /* generic lon/lat grid */

#define GDS_LatLon_nx(gds)	((int) ((gds[6] << 8) + gds[7]))
#define GDS_LatLon_ny(gds)	((int) ((gds[8] << 8) + gds[9]))
#define GDS_LatLon_La1(gds)	INT3(gds[10],gds[11],gds[12])
#define GDS_LatLon_Lo1(gds)	INT3(gds[13],gds[14],gds[15])
#define GDS_LatLon_La2(gds)	INT3(gds[17],gds[18],gds[19])
#define GDS_LatLon_Lo2(gds)	INT3(gds[20],gds[21],gds[22])

#define GDS_LatLon_dx(gds)      INT2(gds[23],gds[24])
#define GDS_LatLon_dy(gds)      INT2(gds[25],gds[26])
#define GDS_Gaussian_nlat(gds)       ((gds[25]<<8)+gds[26])


#define GDS_LatLon_Scan(gds)	(gds[27])

#define GDS_Polar_nx(gds)	((gds[6] << 8) + gds[7])
#define GDS_Polar_ny(gds)	((gds[8] << 8) + gds[9])
#define GDS_Polar_La1(gds)	(INT3(gds[10],gds[11],gds[12]))
#define GDS_Polar_Lo1(gds)	(INT3(gds[13],gds[14],gds[15]))
#define GDS_Polar_Lov(gds)	(INT3(gds[17],gds[18],gds[19]))
#define GDS_Polar_Scan(gds)	(gds[27])
#define GDS_Polar_Dx(gds)	(INT3(gds[20], gds[21], gds[22]))
#define GDS_Polar_Dy(gds)	(INT3(gds[23], gds[24], gds[25]))
#define GDS_Polar_pole(gds)	((gds[26] & 128) == 128)


unsigned char *seek_grib(FILE *file, long *pos, long *len_grib, 
        unsigned char *buffer, unsigned int buf_len);

int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer);

double ibm2flt(unsigned char *ibm);
 
void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
        int n_bits, int n, double ref, double scale);
 
double int_power(double x, int y);

int flt2ieee(float x, unsigned char *ieee);

int wrtieee(float *array, int n, int header, FILE *output);

void levels(int, int);
 
void PDStimes(int time_range, int p1, int p2, int time_unit);

int missing_points(unsigned char *bitmap, int n);

void EC_ext(unsigned char *pds, char *prefix, char *suffix);
/* version 3.1 of grib headers  w. ebisuzaki */
/* this version is incomplete */

#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) * (int) (((a & 0x7f) << 8) + b))
#endif

#define PDS_Len1(pds)		(pds[0])
#define PDS_Len2(pds)		(pds[1])
#define PDS_Len3(pds)		(pds[2])
#define PDS_LEN(pds)		((int) ((pds[0]<<16)+(pds[1]<<8)+pds[2]))
#define PDS_Vsn(pds)		(pds[3])
#define PDS_Center(pds)		(pds[4])
#define PDS_Model(pds)		(pds[5])
#define PDS_Grid(pds)		(pds[6])
#define PDS_HAS_GDS(pds)	((pds[7] & 128) != 0)
#define PDS_HAS_BMS(pds)	((pds[7] & 64) != 0)
#define PDS_PARAM(pds)		(pds[8])
#define PDS_L_TYPE(pds)		(pds[9])
#define PDS_LEVEL1(pds)		(pds[10])
#define PDS_LEVEL2(pds)		(pds[11])

#define PDS_KPDS5(pds)		(pds[8])
#define PDS_KPDS6(pds)		(pds[9])
#define PDS_KPDS7(pds)		((int) ((pds[10]<<8) + pds[11]))

#define PDS_Year(pds)		(pds[12])
#define PDS_Month(pds)		(pds[13])
#define PDS_Day(pds)		(pds[14])
#define PDS_Hour(pds)		(pds[15])
#define PDS_Minute(pds)		(pds[16])
#define PDS_ForecastTimeUnit(pds)	(pds[17])
#define PDS_P1(pds)		(pds[18])
#define PDS_P2(pds)		(pds[19])
#define PDS_TimeRange(pds)	(pds[20])
#define PDS_NumAve(pds)		((int) ((pds[21]<<8)+pds[22]))
#define PDS_NumMissing(pds)	(pds[23])
#define PDS_Century(pds)	(pds[24])
#define PDS_Subcenter(pds)	(pds[25])
#define PDS_DecimalScale(pds)	INT2(pds[26],pds[27])

/* various centers */
#define NMC			7
#define ECMWF			98

/* ECMWF Extensions */

#define PDS_EcLocalId(pds)	(PDS_LEN(pds) >= 41 ? (pds[40]) : 0)
#define PDS_EcClass(pds)	(PDS_LEN(pds) >= 42 ? (pds[41]) : 0)
#define PDS_EcType(pds)		(PDS_LEN(pds) >= 43 ? (pds[42]) : 0)
#define PDS_EcStream(pds)	(PDS_LEN(pds) >= 45 ? (INT2(pds[43], pds[44])) : 0)



#define VERSION "wgrib v1.4 (3-18-96) Wesley Ebisuzaki"

#define CHECK_GRIB

/*
 * wgrib.c extract/inventory grib records
 *
 *                              Wesley Ebisuzaki
 *
 * 11/94 - v1.0
 * 11/94 - v1.1: arbitary size grids, -i option
 * 11/94 - v1.2: bug fixes, ieee option, more info
 * 1/95  - v1.2.4: fix headers for SUN acc
 * 2/95  - v1.2.5: add num_ave in -s listing
 * 2/95  - v1.2.6: change %d to %ld
 * 2/95  - v1.2.7: more output, added some polar stereographic support
 * 2/95  - v1.2.8: max min format changed %f to %g, tidying up more info
 * 3/95  - v1.3.0: fix bug with bitmap, allow numbers > UNDEFINED
 * 3/95  - v1.3.1: print number of missing points (verbose)
 * 3/95  - v1.3.2: -append option added
 * 4/95  - v1.3.2a,b: more output, polar stereo support (-V option)
 * 4/95  - v1.3.3: added ECMWF parameter table (prelim)
 * 6/95  - v1.3.4: nxny from BDS rather than gds?
 * 9/95  - v1.3.4d: speedup in grib write
 * 11/95 - v1.3.4f: new ECMWF parameter table (from Mike Fiorino), EC logic
 * 2/96  - v1.3.4g-h: prelim fix for GDS-less grib files
 * 2/96  - v1.3.4i: faster missing(), -V: "pos n" -> "n" (field 2)
 * 3/96  - v1.4: fix return code (!inventory), and short records near EOF
 *
 */

/*
 * MSEEK = I/O buffer size for seek_grib
 */

#define MSEEK 4096
#define BUFF_ALLOC0	40000


#ifndef min
   #define min(a,b)  ((a) < (b) ? (a) : (b))
   #define max(a,b)  ((a) < (b) ? (b) : (a))
#endif

int main(int argc, char **argv) {

    unsigned char *buffer;
    float *array;
    double temp, rmin, rmax;
    int i, j, nx, ny, file_arg;
    int ndx;
    long int len_grib, pos = 0, nxny, buffer_size, n_dump, count = 1;
    unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;
    FILE *input, *dump_file;
    char line[200];
    enum {BINARY, TEXT, IEEE, GRIB, NONE} output_type = NONE;
    enum {DUMP_RECORD, DUMP_POSITION, DUMP_LIST, INVENTORY} mode = INVENTORY;
    int header = 1;
    long int dump = -1;
    int verbose = 0, append = 0;
    char *dump_file_name = "dump", open_parm[3];
    int return_code = 0;

    if (argc == 1) {
	fprintf(stderr, VERSION);
	fprintf(stderr, "\nPortable Grib decoder for NMC Reanalysis\n");
	fprintf(stderr, "   it slices, dices\n\n");
	fprintf(stderr, " usage: %s [grib file] [options]\n\n", argv[0]);
	fprintf(stderr, "  options:\n");
	fprintf(stderr, "  -v  verbose\n");
	fprintf(stderr, "  -V  even more verbose\n");
	fprintf(stderr, "  -s  short inventory\n");
	fprintf(stderr, "      (default) regular inventory\n\n");

	fprintf(stderr, "  -d [record number]  dump record number\n");
	fprintf(stderr, "  -p [byte position]  dump record at position\n");
	fprintf(stderr, "  -i  dump controlled by stdin (inventory list)\n\n");

	fprintf(stderr, "  -text    dump will be a text file\n");
	fprintf(stderr, "  -ieee    dump will be an ieee file\n");
	fprintf(stderr, "  -grib    dump will be a GRIB file\n");
	fprintf(stderr, "  -bin     dump will be in binary (default)\n");
	fprintf(stderr, "  -nh      dump will have no headers\n");
	fprintf(stderr, "  -h       dump will have headers (default)\n");
	fprintf(stderr, "  -append  to dump file\n");
	fprintf(stderr, "  -o [file] changes name of dump file from 'dump'\n");
	fprintf(stderr, "     will scan if without -d, -p, or -i option\n");
	exit(8);
    }
    file_arg = 0;
    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i],"-v") == 0) {
	    verbose = 1;
	    continue;
	}
	if (strcmp(argv[i],"-V") == 0) {
	    verbose = 2;
	    continue;
	}
	if (strcmp(argv[i],"-s") == 0) {
	    verbose = -1;
	    continue;
	}
	if (strcmp(argv[i],"-text") == 0) {
	    output_type = TEXT;
	    continue;
	}
	if (strcmp(argv[i],"-bin") == 0) {
	    output_type = BINARY;
	    continue;
	}
	if (strcmp(argv[i],"-ieee") == 0) {
	    output_type = IEEE;
	    continue;
	}
	if (strcmp(argv[i],"-grib") == 0) {
	    output_type = GRIB;
	    continue;
	}
	if (strcmp(argv[i],"-nh") == 0) {
	    header = 0;
	    continue;
	}
	if (strcmp(argv[i],"-h") == 0) {
	    header = 1;
	    continue;
	}
	if (strcmp(argv[i],"-append") == 0) {
	    append = 1;
	    continue;
	}
	if (strcmp(argv[i],"-d") == 0) {
	    dump = atol(argv[i+1]);
	    i++;
	    if (output_type == NONE) output_type = BINARY;
	    mode = DUMP_RECORD;
	    continue;
	}
	if (strcmp(argv[i],"-p") == 0) {
	    pos = atol(argv[i+1]);
	    i++;
	    dump = 1;
	    if (output_type == NONE) output_type = BINARY;
	    mode = DUMP_POSITION;
	    continue;
	}
	if (strcmp(argv[i],"-i") == 0) {
	    if (output_type == NONE) output_type = BINARY;
	    mode = DUMP_LIST;
	    continue;
	}
	if (strcmp(argv[i],"-o") == 0) {
	    dump_file_name = argv[i+1];
	    i++;
	    continue;
	}
	if (file_arg == 0) {
	    file_arg = i;
	}
	else {
	    fprintf(stderr,"argument: %s ????\n", argv[i]);
	}
    }


    if ((input = fopen(argv[file_arg],"rb")) == NULL) {
        fprintf(stderr,"could not open file: %s\n", argv[file_arg]);
        exit(7);
    }

    if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL) {
	fprintf(stderr,"not enough memory\n");
    }
    buffer_size = BUFF_ALLOC0;

    /* open output file */
    if (mode != INVENTORY) {
	open_parm[0] = 'w'; open_parm[1] = 'b'; open_parm[2] = '\0';
	if (append) open_parm[0] = 'a';
	if (output_type == TEXT) open_parm[1] = '\0';

	if ((dump_file = fopen(dump_file_name,open_parm)) == NULL) {
	    fprintf(stderr,"could not open dump file\n");
	    exit(8);
        }
    }

    /* skip dump - 1 records */
    for (i = 1; i < dump; i++) {
	msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
	if (msg == NULL) {
	    fprintf(stderr, "ran out of data or bad file\n");
	    exit(8);
	}
	pos += len_grib;
    }
    if (dump > 0) count += dump - 1;
    n_dump = 0;

    for (;;) {
	if (n_dump == 1 && (mode == DUMP_RECORD || mode == DUMP_POSITION)) break;
	if (mode == DUMP_LIST) {
	    if (fgets(line,sizeof(line), stdin) == NULL) break;
            line[sizeof(line) - 1] = 0;
            if (sscanf(line,"%ld:%ld:", &count, &pos) != 2) {
		fprintf(stderr,"bad input from stdin\n");
	        exit(8);
	    }
	}

	msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
	if (msg == NULL) {
	    if (mode == INVENTORY) break;
	    fprintf(stderr,"missing GRIB record(s)\n");
	    exit(8);
	}

        /* read all whole grib record */
        if (len_grib + msg - buffer > buffer_size) {
            buffer_size = len_grib + msg - buffer + 1000;
            buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
            if (buffer == NULL) {
                fprintf(stderr,"ran out of memory\n");
                exit(8);
            }
        }
        read_grib(input, pos, len_grib, buffer);

	/* parse grib message */

	msg = buffer;
        pds = (msg + 8);
        pointer = pds + PDS_LEN(pds);

        if (PDS_HAS_GDS(pds)) {
            gds = pointer;
            pointer += GDS_LEN(gds);
        }
        else {
            gds = NULL;
        }

        if (PDS_HAS_BMS(pds)) {
            bms = pointer;
            pointer += BMS_LEN(bms);
        }
        else {
            bms = NULL;
        }

        bds = pointer;
        pointer += BDS_LEN(bds);

        /* end section - "7777" in ascii */
        if (pointer[0] != 0x37 || pointer[1] != 0x37 ||
            pointer[2] != 0x37 || pointer[3] != 0x37) {
            fprintf(stderr,"\n\n    missing end section\n");
            fprintf(stderr, "%2x %2x %2x %2x\n", pointer[0], pointer[1], 
		pointer[2], pointer[3]);
	    exit(8);
        }

	/* figure out size of array */
	if (gds != NULL) {
	    /* this doesn't work for spherical harmonics */
	    nx = GDS_LatLon_nx(gds);
	    ny = GDS_LatLon_ny(gds);
	}
	else if (bms != NULL) {
	    nx = BMS_nxny(bms);
	    ny = 1;
	}
	else {
	    nx = BDS_NValues(bds);
	    ny = 1;
	}

#ifdef CHECK_GRIB
        if (BDS_NumBits(bds) != 0) {
	    nxny = BDS_NValues(bds);
	    if (bms != NULL) {
		nxny += missing_points(BMS_bitmap(bms),nx*ny);
	    }
	    if (nxny != nx*ny && (gds != NULL && !GDS_Harmonic(gds)) ) {
printf("qqq %d %d\n",GDS_Harmonic(gds));
	        fprintf(stderr,"grib header at record %ld: two values of nxny %ld %d\n",
			count,nxny,nx*ny);
	        fprintf(stderr,"   LEN %d DataStart %d UnusedBits %d #Bits %d\n",
		    BDS_LEN(bds), BDS_DataStart(bds),BDS_UnusedBits(bds),
		    BDS_NumBits(bds));
                return_code = 15;
	    }
	}
#endif
 
        nxny = nx * ny;
	if (verbose <= 0) {
	    printf("%ld:%ld:d=%2.2d%2.2d%2.2d%2.2d:%s:", count, pos, 
                PDS_Year(pds), PDS_Month(pds), PDS_Day(pds), PDS_Hour(pds),
	        k5toa_def(PDS_PARAM(pds),PDS_Center(pds)));
	    if (verbose == 0) printf("kpds5=%d:kpds6=%d:kpds7=%d:TR=%d:P1=%d:P2=%d:TimeU=%d:",
		PDS_PARAM(pds),PDS_KPDS6(pds),PDS_KPDS7(pds),
	        PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                PDS_ForecastTimeUnit(pds));
	    levels(PDS_KPDS6(pds), PDS_KPDS7(pds));
	    printf(":");
	    PDStimes(PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                PDS_ForecastTimeUnit(pds));
	    if (PDS_Center(pds) == ECMWF) EC_ext(pds,"",":");
	    printf("NAve=%d",PDS_NumAve(pds));

	    printf("\n");
	}

        if (verbose > 0) {
	    printf("rec %ld:%ld:date %2.2d%2.2d%2.2d%2.2d %s kpds5=%d "
		"kpds6=%d kpds7=%d levels=(%d,%d) grid=%d ", count, pos, 
                PDS_Year(pds), PDS_Month(pds), PDS_Day(pds), PDS_Hour(pds),
	        k5toa_def(PDS_PARAM(pds),PDS_Center(pds)),PDS_PARAM(pds),
                PDS_KPDS6(pds), PDS_KPDS7(pds), PDS_LEVEL1(pds), 
                PDS_LEVEL2(pds), PDS_Grid(pds));
	    levels(PDS_KPDS6(pds),PDS_KPDS7(pds));
	    printf(" ");
	    if (PDS_Center(pds) == ECMWF) EC_ext(pds,""," ");
	    PDStimes(PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                 PDS_ForecastTimeUnit(pds));
	    if (bms != NULL) 
		printf(" bitmap missing %d", missing_points(BMS_bitmap(bms),nx*ny));
	    printf("\n");

            printf("  timerange %d P1 %d P2 %d  nx %d ny %d GDS grid %d "
		"num_in_ave %d missing %d\n", 
	        PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds), nx, ny, 
                gds == NULL ? -1 : GDS_DataType(gds), 
                PDS_NumAve(pds), PDS_NumMissing(pds));
	}
	if (verbose > 1) {
	    printf("  center %d subcenter %d process %d\n", 
		PDS_Center(pds),PDS_Subcenter(pds),PDS_Model(pds));
	    if (gds && GDS_LatLon(gds)) 
		printf("  latlon: lat  %f to %f by %f  nxny %ld\n"
                       "          long %f to %f by %f, (%d x %d) scan %d"
			" bdsgrid %d\n",
		  0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		  0.001*GDS_LatLon_dy(gds), nxny, 0.001*GDS_LatLon_Lo1(gds),
		  0.001*GDS_LatLon_Lo2(gds), 0.001*GDS_LatLon_dx(gds),
	    	  nx, ny, GDS_LatLon_Scan(gds),
		  BDS_Grid(bds));
/*mf ------- generic lon-lat grid -------- mf*/

	    if (gds && GDS_Generic(gds))  { 
		printf("  GENERIC: lat  %f to %f\n"
                       "           long %f to %f, (%d x %d)",
		       0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		       0.001*GDS_LatLon_Lo1(gds), 0.001*GDS_LatLon_Lo2(gds), 
		       nx, ny);
		printf("\n     lons:\n");
		ndx=23;
		for(i=0;i<nx;i++) {
		  printf("% 10.3f ",0.001*INT3(gds[ndx],gds[ndx+1],gds[ndx+2]));
		  if((i+1)%7==0) printf("\n");
		  ndx+=3;
		}
		if(nx%7!=0) printf("\n");

		printf("     lats:\n");
		for(j=0;j<ny;j++) {
		  printf("% 10.3f ",0.001*INT3(gds[ndx],gds[ndx+1],gds[ndx+2]));
		  if((j+1)%7==0) printf("\n");
		  ndx+=3;
		}
		if(ny%7!=0) printf("\n");

	      }
	    if (gds && GDS_Gaussian(gds)) 
		printf("  gaussian: lat  %f to %f\n"
                       "            long %f to %f by %f, (%d x %d) scan %d"
			" bdsgrid %d\n",
		  0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		  0.001*GDS_LatLon_Lo1(gds), 0.001*GDS_LatLon_Lo2(gds), 
		  0.001*GDS_LatLon_dx(gds),
	    	  nx, ny, GDS_LatLon_Scan(gds),
		  BDS_Grid(bds));
	    if (gds && GDS_Polar(gds))
		printf("  polar stereo: Lat1 %f Long1 %f Orient %f\n"
			"     %s pole (%d x %d) Dx %d Dy %d scan %d\n",
		    0.001*GDS_Polar_La1(gds),0.001*GDS_Polar_Lo1(gds),
		    0.001*GDS_Polar_Lov(gds),
		    GDS_Polar_pole(gds) == 0 ? "north" : "south", nx,ny,
		    GDS_Polar_Dx(gds),GDS_Polar_Dy(gds),
		    GDS_Polar_Scan(gds));

	    printf("\n");
	}

	if (mode != INVENTORY && output_type == GRIB) {
	    fwrite((void *) msg, sizeof(char), len_grib, dump_file);
	    n_dump++;
	}

	if ((mode != INVENTORY && output_type != GRIB) || verbose > 1) {
	    /* decode numeric data */
 
            if ((array = (float *) malloc(sizeof(float) * nxny)) == NULL) {
                fprintf(stderr,"memory problems\n");
                exit(8);
            }

	    temp = int_power(10.0, - PDS_DecimalScale(pds));

	    BDS_unpack(array, bds + 11, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
	        temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));

	    if (verbose > 1) {
		rmin = FLT_MAX;
		rmax = -FLT_MAX;
	        for (i = 0; i < nxny; i++) {
		    if (fabs(array[i]-UNDEFINED) > 0.0001*UNDEFINED) {
	                rmin = min(rmin,array[i]);
	                rmax = max(rmax,array[i]);
		    }
	        }
	        printf("  min/max data %g %g  num bits %d "
			" BDS_Ref %g  DecScale %d BinScale %d\n", 
		    rmin, rmax, BDS_NumBits(bds), BDS_RefValue(bds),
		    PDS_DecimalScale(pds), BDS_BinScale(bds));
	    }

	    if (mode != INVENTORY && output_type != GRIB) {
	        if (output_type == BINARY) {
	            i = nxny * sizeof(float);
	            if (header) fwrite((void *) &i, sizeof(int), 1, dump_file);
	            fwrite((void *) array, sizeof(float), nxny, dump_file);
	            if (header) fwrite((void *) &i, sizeof(int), 1, dump_file);
	        }
		else if (output_type == IEEE) {
		    wrtieee(array, nxny, header, dump_file);
		}
	        else if (output_type == TEXT) {
	            /* number of points in grid */
	            if (header) fprintf(dump_file, "%d %d\n", nx, ny);

	            for (i = 0; i < nxny; i++) {
		        fprintf(dump_file,"%g\n", array[i]);
		    }
	        }
	        n_dump++;
	    }
	    free(array);
	    if (verbose > 0) printf("\n");
	}
	    
        pos += len_grib;
        count++;
    }

    if (mode != INVENTORY) {
	if (ferror(dump_file)) {
		fprintf(stderr,"error writing %s\n",dump_file_name);
		exit(8);
	}
    }
    fclose(input);
    return (return_code);
}
/*
 * find next grib header
 *
 * file = what do you think?
 * pos = initial position to start looking at  ( = 0 for 1st call)
 *       returns with position of next grib header (units=bytes)
 * len_grib = length of the grib record (bytes)
 * buffer[buf_len] = buffer for reading/writing
 *
 * returns (char *) to start of GRIB header+PDS
 *         NULL if not found
 *
 * adapted from SKGB (Mark Iredell)
 *
 * v1.1 9/94 Wesley Ebisuzaki
 * v1.2 3/96 Wesley Ebisuzaki handles short records at end of file
 *
 */

#ifndef min
   #define min(a,b)  ((a) < (b) ? (a) : (b))
#endif

#define NTRY 3
/* #define LEN_HEADER_PDS (28+42+100) */
#define LEN_HEADER_PDS (28+8)

unsigned char *seek_grib(FILE *file, long *pos, long *len_grib, 
        unsigned char *buffer, unsigned int buf_len) {

    int i, j, len;

    for (j = 0; j < NTRY; j++) {

        if (fseek(file, *pos, SEEK_SET) == -1) {
            *len_grib = 0;
            return (unsigned char *) NULL;
        }
     
        i = fread(buffer, sizeof (unsigned char), buf_len, file);
     
        len = i - LEN_HEADER_PDS;
     
        for (i = 0; i < len; i++) {
            if (buffer[i] == 'G' && buffer[i+1] == 'R' && buffer[i+2] == 'I'
                && buffer[i+3] == 'B' && buffer[i+7] == 1) {
                    *pos = i + *pos;
                    *len_grib = (buffer[i+4] << 16) + (buffer[i+5] << 8) +
                            buffer[i+6];
                    return (buffer+i);
            }
        }
	*pos = *pos + (buf_len - LEN_HEADER_PDS);
    }

    *len_grib = 0;
    return (unsigned char *) NULL;
}
/* wesley ebisuzaki v 1.1 */


double ibm2flt(unsigned char *ibm) {

	int positive, power;
	unsigned int abspower;
	long int mant;
	double value, exp;

	positive = (ibm[0] & 0x80) == 0;
	mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
	power = (int) (ibm[0] & 0x7f) - 64;
	abspower = power > 0 ? power : -power;


	/* calc exp */
	exp = 16.0;
	value = 1.0;
	while (abspower) {
		if (abspower & 1) {
			value *= exp;
		}
		exp = exp * exp;
		abspower >>= 1;
	}

	if (power < 0) value = 1.0 / value;
	value = value * mant / 16777216.0;
	if (positive == 0) value = -value;
	return value;
}
	
/*
 * read_grib.c
 *
 * v1.0 9/94 Wesley Ebisuzaki
 *
 */

int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer) {

    int i;


    if (fseek(file, pos, SEEK_SET) == -1) {
	    return 0;
    }

    i = fread(buffer, sizeof (unsigned char), len_grib, file);
    return (i == len_grib);
}

/*
 * w. ebisuzaki
 *
 *  return x**y
 *
 *
 *  input: double x
 *	   int y
 */
double int_power(double x, int y) {

	double value;

	if (y < 0) {
		y = -y;
		x = 1.0 / x;
	}
	value = 1.0;

	while (y) {
		if (y & 1) {
			value *= x;
		}
		x = x * x;
		y >>= 1;
	}
	return value;
}

/* Wesley Ebisuzaki */

static char *cname5[256] = {
  "var0",  /* 0 */
  "PRES",  /* 1 */
  "PRMSL",  /* 2 */
  "PTEND",  /* 3 */
  "var4",  /* 4 */
  "var5",  /* 5 */
  "GP",  /* 6 */
  "HGT",  /* 7 */
  "DIST",  /* 8 */
  "HSTDV",  /* 9 */
  "HVAR",  /* 10 */
  "TMP",  /* 11 */
  "VTMP",  /* 12 */
  "POT",  /* 13 */
  "EPOT",  /* 14 */
  "TMAX",  /* 15 */
  "TMIN",  /* 16 */
  "DPT",  /* 17 */
  "DEPR",  /* 18 */
  "LAPR",  /* 19 */
  "VISIB",  /* 20 */
  "RDSP1",  /* 21 */
  "RDSP2",  /* 22 */
  "RDSP3",  /* 23 */
  "var24",  /* 24 */
  "TMPA",  /* 25 */
  "PRESA",  /* 26 */
  "GPA",  /* 27 */
  "WVSP1",  /* 28 */
  "WVSP2",  /* 29 */
  "WVSP3",  /* 30 */
  "WDIR",  /* 31 */
  "WIND",  /* 32 */
  "UGRD",  /* 33 */
  "VGRD",  /* 34 */
  "STRM",  /* 35 */
  "VPOT",  /* 36 */
  "MNTSF",  /* 37 */
  "SGCVV",  /* 38 */
  "VVEL",  /* 39 */
  "DZDT",  /* 40 */
  "ABSV",  /* 41 */
  "ABSD",  /* 42 */
  "RELV",  /* 43 */
  "RELD",  /* 44 */
  "VUCSH",  /* 45 */
  "VVCSH",  /* 46 */
  "DIRC",  /* 47 */
  "SPC",  /* 48 */
  "UOGRD",  /* 49 */
  "VOGRD",  /* 50 */
  "SPFH",  /* 51 */
  "RH",  /* 52 */
  "MIXR",  /* 53 */
  "PWAT",  /* 54 */
  "VAPP",  /* 55 */
  "SATD",  /* 56 */
  "EVP",  /* 57 */
  "CICE",  /* 58 */
  "PRATE",  /* 59 */
  "TSTM",  /* 60 */
  "APCP",  /* 61 */
  "NCPCP",  /* 62 */
  "ACPCP",  /* 63 */
  "SRWEQ",  /* 64 */
  "WEASD",  /* 65 */
  "SNOD",  /* 66 */
  "MIXHT",  /* 67 */
  "TTHDP",  /* 68 */
  "MTHD",  /* 69 */
  "MTHA",  /* 70 */
  "TCDC",  /* 71 */
  "CDCON",  /* 72 */
  "LCDC",  /* 73 */
  "MCDC",  /* 74 */
  "HCDC",  /* 75 */
  "CWAT",  /* 76 */
  "var77",  /* 77 */
  "SNOC",  /* 78 */
  "SNOL",  /* 79 */
  "WTMP",  /* 80 */
  "LAND",  /* 81 */
  "DSLM",  /* 82 */
  "SFCR",  /* 83 */
  "ALBDO",  /* 84 */
  "TSOIL",  /* 85 */
  "SOILM",  /* 86 */
  "VEG",  /* 87 */
  "SALTY",  /* 88 */
  "DEN",  /* 89 */
  "RUNOF",  /* 90 */
  "ICEC",  /* 91 */
  "ICETK",  /* 92 */
  "DICED",  /* 93 */
  "SICED",  /* 94 */
  "UICE",  /* 95 */
  "VICE",  /* 96 */
  "ICEG",  /* 97 */
  "ICED",  /* 98 */
  "SNOM",  /* 99 */
  "HTSGW",  /* 100 */
  "WVDIR",  /* 101 */
  "WVHGT",  /* 102 */
  "WVPER",  /* 103 */
  "SWDIR",  /* 104 */
  "SWELL",  /* 105 */
  "SWPER",  /* 106 */
  "DIRPW",  /* 107 */
  "PERPW",  /* 108 */
  "DIRSW",  /* 109 */
  "PERSW",  /* 110 */
  "NSWRS",  /* 111 */
  "NLWRS",  /* 112 */
  "NSWRT",  /* 113 */
  "NLWRT",  /* 114 */
  "LWAVR",  /* 115 */
  "SWAVR",  /* 116 */
  "GRAD",  /* 117 */
  "var118",  /* 118 */
  "var119",  /* 119 */
  "var120",  /* 120 */
  "LHTFL",  /* 121 */
  "SHTFL",  /* 122 */
  "BLYDP",  /* 123 */
  "UFLX",  /* 124 */
  "VFLX",  /* 125 */
  "WMIXE",  /* 126 */
  "IMGD",  /* 127 */
  "MSLSA",  /* 128 */
  "MSLMA",  /* 129 */
  "MSLET",  /* 130 */
  "LFTX",  /* 131 */
  "4LFTX",  /* 132 */
  "KX",  /* 133 */
  "SX",  /* 134 */
  "MCONV",  /* 135 */
  "VSSH",  /* 136 */
  "TSLSA",  /* 137 */
  "BVF2",  /* 138 */
  "PVMW",  /* 139 */
  "CRAIN",  /* 140 */
  "CFRZR",  /* 141 */
  "CICEP",  /* 142 */
  "CSNOW",  /* 143 */
  "SOILW",  /* 144 */
  "PEVPR",  /* 145 */
  "CWORK",  /* 146 */
  "U-GWD",  /* 147 */
  "V-GWD",  /* 148 */
  "PV___",  /* 149 */
  "var150",  /* 150 */
  "var151",  /* 151 */
  "var152",  /* 152 */
  "MFXDV",  /* 153 */
  "var154",  /* 154 */
  "GFLUX",  /* 155 */
  "CIN",  /* 156 */
  "CAPE",  /* 157 */
  "TKE",  /* 158 */
  "CONDP",  /* 159 */
  "CSUSF",  /* 160 */
  "CSDSF",  /* 161 */
  "CSULF",  /* 162 */
  "CSDLF",  /* 163 */
  "CFNSF",  /* 164 */
  "CFNLF",  /* 165 */
  "VBDSF",  /* 166 */
  "VDDSF",  /* 167 */
  "NBDSF",  /* 168 */
  "NDDSF",  /* 169 */
  "USTR",  /* 170 */
  "VSTR",  /* 171 */
  "MFLX",  /* 172 */
  "LMH",  /* 173 */
  "LMV",  /* 174 */
  "SGLYR",  /* 175 */
  "NLAT",  /* 176 */
  "NLON",  /* 177 */
  "UMAS",  /* 178 */
  "VMAS",  /* 179 */
  "var180",  /* 180 */
  "LPSX",  /* 181 */
  "LPSY",  /* 182 */
  "HGTX",  /* 183 */
  "HGTY",  /* 184 */
  "STDZ",  /* 185 */
  "STDU",  /* 186 */
  "STDV",  /* 187 */
  "STDQ",  /* 188 */
  "STDT",  /* 189 */
  "CBUW",  /* 190 */
  "CBVW",  /* 191 */
  "CBUQ",  /* 192 */
  "CBVQ",  /* 193 */
  "CBTW",  /* 194 */
  "CBQW",  /* 195 */
  "CBMZW",  /* 196 */
  "CBTZW",  /* 197 */
  "CBTMW",  /* 198 */
  "STDRH",  /* 199 */
  "SDTZ",  /* 200 */
  "ICWAT",  /* 201 */
  "SDTU",  /* 202 */
  "SDTV",  /* 203 */
  "DSWRF",  /* 204 */
  "DLWRF",  /* 205 */
  "SDTQ",  /* 206 */
  "MSTAV",  /* 207 */
  "SFEXC",  /* 208 */
  "MIXLY",  /* 209 */
  "SDTT",  /* 210 */
  "USWRF",  /* 211 */
  "ULWRF",  /* 212 */
  "CDLYR",  /* 213 */
  "CPRAT",  /* 214 */
  "TTDIA",  /* 215 */
  "TTRAD",  /* 216 */
  "TTPHY",  /* 217 */
  "PREIX",  /* 218 */
  "TSD1D",  /* 219 */
  "NLSGP",  /* 220 */
  "SDTRH",  /* 221 */
  "5WAVH",  /* 222 */
  "CWAT",  /* 223 */
  "PLTRS",  /* 224 */
  "var225",  /* 225 */
  "BMIXL",  /* 226 */
  "AMIXL",  /* 227 */
  "PEVAP",  /* 228 */
  "SNOHF",  /* 229 */
  "var230",  /* 230 */
  "MFLUX",  /* 231 */
  "DTRF",  /* 232 */
  "UTRF",  /* 233 */
  "BGRUN",  /* 234 */
  "SSRUN",  /* 235 */
  "var236",  /* 236 */
  "OZONE",  /* 237 */
  "SNOC",  /* 238 */
  "SNOT",  /* 239 */
  "GLCR",  /* 240 */
  "LRGHR",  /* 241 */
  "CNVHR",  /* 242 */
  "CNVMR",  /* 243 */
  "SHAHR",  /* 244 */
  "SHAMR",  /* 245 */
  "VDFHR",  /* 246 */
  "VDFUA",  /* 247 */
  "VDFVA",  /* 248 */
  "VDFMR",  /* 249 */
  "SWHR",  /* 250 */
  "LWHR",  /* 251 */
  "CD",  /* 252 */
  "FRICV",  /* 253 */
  "RI",  /* 254 */
  "var255",  /* 255 */
};

extern char *cname5_ecmwf[256];

char *k5toa_def(int i, int center) {

    if (center == ECMWF) return cname5_ecmwf[i];
    return cname5[i];
}

int atok5_def(char *name, int center) {

    int i;

    if (center == ECMWF) {
        for (i = 0; i <= 255; i++) {
            if (strcmp(name, cname5_ecmwf[i]) == 0) return i;
        }
        return -1;
    }

    for (i = 0; i <= 255; i++) {
        if (strcmp(name,cname5[i]) == 0) return i;
    }
    return -1;
}

/* for simple unpacking of a grid */
/* wesley ebisuzaki */

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
	int n_bits, int n, double ref, double scale) {

    int i, j, k;
    unsigned int map_mask, bit_mask;

    map_mask = bit_mask = 128;

    for (i = 0; i < n; i++) {
	if (bitmap) {
	    j = (*bitmap & map_mask);
	    if ((map_mask >>= 1) == 0) {
		map_mask = 128;
		bitmap++;
	    }
	    if (j == 0) {
		*flt++ = UNDEFINED;
		continue;
	    }
	}

	j = 0;
	k = n_bits;
	while (k) {
	    if (k >= 8 && bit_mask == 128) {
		j = 256 * j + *bits;
		bits++;
		k -= 8;
	    }
	    else {
	        j = j + j + ((*bits & bit_mask) != 0);
		if ((bit_mask >>= 1) == 0) {
		    bits++;
		    bit_mask = 128;
		}
		k--;
	    }
	}
	*flt++ = ref + scale*j;
   }
   return;
}

/*
 * convert a float to an ieee single precision number v1.0
 * (big endian)
 *                      Wesley Ebisuzaki
 *
 * bugs: doesn't handle subnormal numbers
 */

int flt2ieee(float x, unsigned char *ieee) {

	int sign, exp, i;
	double mant;

	if (x == 0.0) {
		ieee[0] = ieee[1] = ieee[2] = ieee[3] = 0;
		return 0;
	}

	/* sign bit */
	if (x < 0.0) {
		sign = 128;
		x = -x;
	}
	else sign = 0;
	mant = frexp((double) x, &exp);

	/* round up by adding 2**-24 */
	mant = mant + 1.0/16777216.0;
	if (mant >= 1.0) {
		mant = 0.5;
		exp++;
	}

	mant = mant + mant - 1.0;
	exp = exp - 1 + 127;

	if (exp < 0) {
		/* signed zero */
		ieee[0] = sign;
		ieee[1] = ieee[2] = ieee[3] = 0;
		return 0;
	}
	if (exp > 255) {
		/* signed infinity */
		ieee[0] = sign + 127;
		ieee[1] = 128;
                ieee[2] = ieee[3] = 0;
                return 0;
	}
	/* normal number */

	ieee[0] = sign + (exp >> 1);

	mant = mant * 128.0;
	i = floor(mant);
	mant = mant - i;
	ieee[1] = ((exp & 1) << 7) | i;

	mant = mant * 256.0;
	i = floor(mant);
	mant = mant - i;
	ieee[2] = i;

	ieee[3] = floor(mant*256.0);

	return 0;
}


/* wesley ebisuzaki v1.0
 *
 * write ieee file -- big endin format
 *
 * input float *array		data to be written
 *	 int n			size of array
 *	 int header		1 for f77 style header 0 for none
 *				(header is 4 byte header
 *	 FILE *output		output file
 */

int wrtieee(float *array, int n, int header, FILE *output) {

	unsigned long int l;
	int i;
	unsigned char c4[4], h4[4];

	if (header) {
		l = n * 4;
		for (i = 0; i < 4; i++) {
			h4[i] = l & 255;
			l >>= 8;
		}
		putc(h4[3],output);
		putc(h4[2],output);
		putc(h4[1],output);
		putc(h4[0],output);

	}
	for (i = 0; i < n; i++) {
		flt2ieee(array[i], c4);
		putc(c4[0],output);
		putc(c4[1],output);
		putc(c4[2],output);
		putc(c4[3],output);
	}
	if (header) {
		putc(h4[3],output);
		putc(h4[2],output);
		putc(h4[1],output);
		putc(h4[0],output);
	}
	return 0;
}


/* wesley ebisuzaki v1.0
 *
 * levels.c
 *
 * prints out a simple description of kpds6, kpds7
 *    (level/layer data)
 *  kpds6 = octet 10 of the PDS
 *  kpds7 = octet 11 and 12 of the PDS
 *    (kpds values are from NMC's grib routines)
 *
 * the description of the levels is 
 *   (1) incomplete
 *   (2) include some NMC-only values (>= 200?)
 */

void levels(int kpds6, int kpds7) {

	int o11, o12;

	/* octets 11 and 12 */
	o11 = kpds7 / 256;
	o12 = kpds7 % 256;


	switch (kpds6) {

	case 1: printf("sfc");
		break;
	case 2: printf("cld base");
		break;
	case 3: printf("cld top");
		break;
	case 4: printf("0C isotherm");
		break;
	case 5: printf("cond lev");
		break;
	case 6: printf("max wind lev");
		break;
	case 7: printf("tropopause");
		break;
	case 8: printf("nom. top");
		break;
	case 9: printf("sea bottom");
		break;
	case 200:
	case 10: printf("atmos col");
		break;

	case 12:
	case 212: printf("low cld bot");
		break;
	case 13:
	case 213: printf("low cld top");
		break;
	case 14:
	case 214: printf("low cld lay");
		break;
	case 22:
	case 222: printf("mid cld bot");
		break;
	case 23:
	case 223: printf("mid cld top");
		break;
	case 24:
	case 224: printf("mid cld lay");
		break;
	case 32:
	case 232: printf("high cld bot");
		break;
	case 33:
	case 233: printf("high cld top");
		break;
	case 34:
	case 234: printf("high cld lay");
		break;


	case 100: printf("%d mb",kpds7);
	 	break;
	case 101: printf("%d-%d mb",o11*10,o12*10);
	 	break;
	case 102: printf("MSL");
	 	break;
	case 103: printf("%d m above MSL",kpds7);
	 	break;
	case 104: printf("%d-%d m above msl",o11*100,o12*100);
	 	break;
	case 105: printf("%d m above gnd",kpds7);
	 	break;
	case 106: printf("%d-%d m above gnd",o11*100,o12*100);
	 	break;
	case 107: printf("sigma=%.4f",kpds7/10000.0);
	 	break;
	case 108: printf("sigma %.2f-%.2f",o11/100.0,o12/100.0);
	 	break;
	case 109: printf("hybrid lev %d",kpds7);
	 	break;
	case 110: printf("hybrid %d-%d",o11,o12);
	 	break;
	case 111: printf("%d cm down",kpds7);
	 	break;
	case 112: printf("%d-%d cm down",o11,o12);
	 	break;
	case 113: printf("%dK",kpds7);
	 	break;
	case 114: printf("%d-%dK",475-o11,475-o12);
	 	break;
	case 115: printf("%d mb above gnd",kpds7);
	 	break;
	case 116: printf("%d-%d mb above gnd",o11,o12);
	 	break;
	case 121: printf("%d-%d mb",1100-o11,1100-o12);
	 	break;
	default:
	 	break;
	}
}

/*
 * PDStimes.c   v1.1a wesley ebisuzaki
 *
 * prints something readable for time code in grib file
 *
 * not all cases decoded
 * for NMC/NCAR Reanalysis
 */

static char *units[] = {
	"min", "hr", "d", "mon", "yr",
	"decade", "normal", "century"};


void PDStimes(int time_range, int p1, int p2, int time_unit) {

	char *unit;
	enum {anal, fcst, unknown} type;
	int fcst_len;

	if (time_unit >= 0 && time_unit <= 7) unit = units[time_unit];
	else unit = "";

	/* figure out if analysis or forecast */
	/* in GRIB, there is a difference between init and uninit analyses */
	/* not case at NMC .. no longer run initialization */
	/* ignore diff between init an uninit analyses */

	switch (time_range) {

	case 0:
	case 1:
	case 113:
	case 114:
	case 118:
		if (p1 == 0) type = anal;
		else {
			type = fcst;
			fcst_len = p1;
		}
		break;
	case 10: type = fcst; /* way NMC uses it, should be unknown? */
		fcst_len = p1*256 + p2;
		break;

	case 51:
	case 123:
	case 124:
		type = anal;
		break;

	default: type = unknown;
		break;
	}

	/* ----------------------------------------------- */

	if (type == anal) printf("anl:");
	else if (type == fcst) printf("%d%s fcst:",fcst_len,unit);


	if (time_range == 123 || time_range == 124) {
		if (p1 != 0) printf("start@ %d%s:",p1,unit);
	}


	/* print time range */


	switch (time_range) {

	case 0:
	case 1:
	case 10:
		break;
	case 2: printf("valid %d-%d%s:",p1,p2,unit);
		break;
	case 3: printf("%d-%d%s ave:",p1,p2,unit);
		break;
	case 4: printf("%d-%d%s acc:",p1,p2,unit);
		break;
	case 5: printf("%d-%d%s diff:",p1,p2,unit);
		break;
	case 51: printf("clim mean %d %d%s:",p1,p2,unit);
		break;
	case 113:
	case 123:
		printf("ave@%d%s:",p2,unit);
		break;
	case 114:
	case 124:
		printf("acc@%d%s:",p2,unit);
		break;
	case 115:
		printf("ave of fcst:%d to %d%s:",p1,p2,unit);
		break;
	case 116:
		printf("acc of fcst:%d to %d%s:",p1,p2,unit);
		break;
	case 118: 
		printf("var@%d%s:",p2,unit);
		break;
	default: printf("time?:");
	}
}


/*
 * v1.1   number of missing data points w. ebisuzaki
 *  change from 1.0: just faster my dear
 *
 */

static int mask[] = {128, 64, 32, 16, 8, 4, 2, 1};

int missing_points(unsigned char *bitmap, int n) {

    int i, count;
    unsigned int tmp;
    if (bitmap == NULL) return 0;

    count = 0;
    while (n >= 8) {
	tmp = *bitmap++;
	n -= 8;
	if ((mask[0] & tmp) == 0) count++;
	if ((mask[1] & tmp) == 0) count++;
	if ((mask[2] & tmp) == 0) count++;
	if ((mask[3] & tmp) == 0) count++;
	if ((mask[4] & tmp) == 0) count++;
	if ((mask[5] & tmp) == 0) count++;
	if ((mask[6] & tmp) == 0) count++;
	if ((mask[7] & tmp) == 0) count++;
    }

    for (i = 0; i < n; i++) {
	count += ((mask[i] & *bitmap) == 0);
    }
    return count;
}
/*
 *  parameter table for ECMWF
 *
 *  this table is from Mike Fiorino (LLNL) and was used in
 *  the ECMWF's reanalysis project.  This table was a superset
 *  of the table that I obtained from NCAR.
 *
 */

char *cname5_ecmwf[256] = {
  "var0",  /* 0 */
  "PRES",  /* 1 */
  "PRMSL",  /* 2 */
  "PTEND",  /* 3 */
  "var4",  /* 4 */
  "var5",  /* 5 */
  "GP",  /* 6 */
  "HGT",  /* 7 */
  "DIST",  /* 8 */
  "HSTDV",  /* 9 */
  "HVAR",  /* 10 */
  "TMP",  /* 11 */
  "VTMP",  /* 12 */
  "POT",  /* 13 */
  "EPOT",  /* 14 */
  "TMAX",  /* 15 */
  "TMIN",  /* 16 */
  "DPT",  /* 17 */
  "DEPR",  /* 18 */
  "LAPR",  /* 19 */
  "VISIB",  /* 20 */
  "RDSP1",  /* 21 */
  "RDSP2",  /* 22 */
  "RDSP3",  /* 23 */
  "var24",  /* 24 */
  "TMPA",  /* 25 */
  "PRESA",  /* 26 */
  "GPA",  /* 27 */
  "WVSP1",  /* 28 */
  "WVSP2",  /* 29 */
  "WVSP3",  /* 30 */
  "WDIR",  /* 31 */
  "WIND",  /* 32 */
  "UGRD",  /* 33 */
  "VGRD",  /* 34 */
  "STRM",  /* 35 */
  "VPOT",  /* 36 */
  "MNTSF",  /* 37 */
  "SGCVV",  /* 38 */
  "VVEL",  /* 39 */
  "DZDT",  /* 40 */
  "ABSV",  /* 41 */
  "ABSD",  /* 42 */
  "RELV",  /* 43 */
  "RELD",  /* 44 */
  "VUCSH",  /* 45 */
  "VVCSH",  /* 46 */
  "DIRC",  /* 47 */
  "SPC",  /* 48 */
  "UOGRD",  /* 49 */
  "VOGRD",  /* 50 */
  "SPFH",  /* 51 */
  "RH",  /* 52 */
  "MIXR",  /* 53 */
  "PWAT",  /* 54 */
  "VAPP",  /* 55 */
  "SATD",  /* 56 */
  "EVP",  /* 57 */
  "CICE",  /* 58 */
  "PRATE",  /* 59 */
  "TSTM",  /* 60 */
  "APCP",  /* 61 */
  "NCPCP",  /* 62 */
  "ACPCP",  /* 63 */
  "SRWEQ",  /* 64 */
  "WEASD",  /* 65 */
  "SNOD",  /* 66 */
  "MIXHT",  /* 67 */
  "TTHDP",  /* 68 */
  "MTHD",  /* 69 */
  "MTHA",  /* 70 */
  "TCDC",  /* 71 */
  "CDCON",  /* 72 */
  "LCDC",  /* 73 */
  "MCDC",  /* 74 */
  "HCDC",  /* 75 */
  "CWAT",  /* 76 */
  "var77",  /* 77 */
  "SNOC",  /* 78 */
  "SNOL",  /* 79 */
  "WTMP",  /* 80 */
  "LAND",  /* 81 */
  "DSLM",  /* 82 */
  "SFCR",  /* 83 */
  "ALBDO",  /* 84 */
  "TSOIL",  /* 85 */
  "SOILM",  /* 86 */
  "VEG",  /* 87 */
  "SALTY",  /* 88 */
  "DEN",  /* 89 */
  "RUNOF",  /* 90 */
  "ICEC",  /* 91 */
  "ICETK",  /* 92 */
  "DICED",  /* 93 */
  "SICED",  /* 94 */
  "UICE",  /* 95 */
  "VICE",  /* 96 */
  "ICEG",  /* 97 */
  "ICED",  /* 98 */
  "SNOM",  /* 99 */
  "HTSGW",  /* 100 */
  "WVDIR",  /* 101 */
  "WVHGT",  /* 102 */
  "WVPER",  /* 103 */
  "SWDIR",  /* 104 */
  "SWELL",  /* 105 */
  "SWPER",  /* 106 */
  "DIRPW",  /* 107 */
  "PERPW",  /* 108 */
  "DIRSW",  /* 109 */
  "PERSW",  /* 110 */
  "NSWRS",  /* 111 */
  "NLWRS",  /* 112 */
  "NSWRT",  /* 113 */
  "NLWRT",  /* 114 */
  "LWAVR",  /* 115 */
  "SWAVR",  /* 116 */
  "GRAD",  /* 117 */
  "var118",  /* 118 */
  "var119",  /* 119 */
  "var120",  /* 120 */
  "LHTFL",  /* 121 */
  "SHTFL",  /* 122 */
  "BLYDP",  /* 123 */
  "UFLX",  /* 124 */
  "VFLX",  /* 125 */
  "WMIXE",  /* 126 */

/*
          ECMWF Version of Table 2 for FM92-VIII Ext. GRIB

 (This table is used instead of WMO version of Table 2 of FM92-VIII Ext. GRIB)

CODE FIGURE                 FIELD NAME                        UNITS


   * denotes field accumulated since start of forecast

  abbrev. Mike Fiorino
  note: overlap with standard & NMC abbrev
  note: ECMWF should not have redefined #127

*/

   "at",      /*   127  Atmospheric Tide                     -      */
   "bdv",     /*   128  Budget Values                        -      */
   "zg",      /*   129  Geopotential [m2 s-1]                       */
   "ta",      /*   130  Temperature [K]                             */
   "ua",      /*   131  U-component of Wind [ms-1]                  */
   "va",      /*   132  V-component of Wind [ms-1]                  */
   "hus",     /*   133  Specific Humidity [kg/kg]                   */
   "pss",     /*   134  Surface Pressure [Pa]                       */
   "wa",      /*   135  Vertical Velocity [Pa s-1]                  */
   "prwa",    /*   136  preciptable water (vapor+drops+ice) [m]     */
   "prw",     /*   137  pecipitable water [m]                       */
   "rvort",   /*   138  vorticity [s-1]                             */
   "tso1",    /*   139  soil moisture level 1 [m (H20)]             */
   "mrso1",   /*   140  soil temperature level 1 [K]                */
   "snd",     /*   141  Snow Depth [m]                              */
   "prl",     /*   142  Large Scale Precipitation [m]               */
   "prc",     /*   143  Convective Precipitation [m]                */
   "prs",     /*   144  Snow Fall                                   */
   "bld",     /*   145  Boundary Layer Dissipation [Wm-2]           */
   "hfss",    /*   146  Surface Flux of Sensible Heat [Wm-2]        */
   "hfls",    /*   147  Surface Flux of Latent Heat [Wm-2]          */
   "v148",    /*   148                                              */
   "v149",    /*   149                                              */
   "v150",    /*   150                                              */
   "psl",     /*   151  Mean Sea Level (MSL) Pressure [Pa]    Pa    */
   "logpsl",  /*   152  Log Surface Pressure                  -     */
   "v153",    /*   153                                              */
   "v154",    /*   154                                              */
   "div",     /*   155  Divergence [s-1]                            */
   "zg",      /*   156  Height (Geopotential) [m]                   */
   "hur",     /*   157  Relative Humidity [%]                       */
   "pstn",    /*   158  Tendency of Surface Pressure [Pa s-1]       */
   "v159",    /*   159                                              */
   "v160",    /*   150                                              */
   "v161",    /*   161                                              */
   "v162",    /*   162                                              */
   "v163",    /*   163                                              */
   "clt",     /*   164  cloud cover total [0-1]                     */
   "uas",     /*   165  U-wind at 10 [ms-1]                         */
   "vas",     /*   166  V-wind at 10 [ms-1]                         */
   "tas",     /*   167  Temperature at 2 m [K]                      */
   "tds",     /*   168  Dewpoint at 2 m [K]                         */
   "rsds",    /*   169  Downward SW (sfc) [Wm-2]                    */
   "tso2",    /*   170  soil temperature [K]                        */
   "mrso2",   /*   171  soil wetness level 2 [m (H20)]              */
   "lsm",     /*   172  land-sea mask [(0,1)]                       */
   "sfr",     /*   173  sfc roughness [m]                           */
   "albs",    /*   174  Albedo [0-1]                                */
   "rlds",    /*   175  Downard LW (sfc) [Wm-2]                     */
   "rss",     /*   176  Net Shortwave Radiation (surface) [Wm-2]    */
   "rls",     /*   177  Net Longwave Radiation (surface) [Wm-2]     */
   "rst",     /*   178  Net Shortwave Radiation (toa) [Wm-2         */
   "rlt",     /*   179  Net Longwave Radiation (toa) [Wm-2]         */
   "tauu",    /*   180  U-component of Surface Wind Stress [Nm-2]   */
   "tauv",    /*   181  V-component of Surface Wind Stress [Nm-2]   */
   "evs",     /*   182  Evaporation [m (H2O)]                       */
   "tso3",    /*   183  soil temp level 3 [m (H2O)]                 */
   "mrso3",   /*   184  soil moisture level 3 [K]                   */
   "clcc",    /*   185  cloud convective [0-1]                      */
   "cll",     /*   186  cloud low [0-1]                             */
   "clm",     /*   187  cloud mid [0-1]                             */
   "clh",     /*   188  cloud high [0-1]                            */
   "v189",    /*   189                                              */
   "orgv",    /*   190  orographic variance                         */
   "orgvew",  /*   191  orographic variance e-w                     */
   "orgvns",  /*   192  orographic variance n-s                     */
   "orgvnwse",/*   193  orographic variance nw-se                   */
   "orgvnesw",/*   194  orographic variance ne-sw                   */
   "gwsv",    /*   195  gravity wave stress n-s                     */
   "gwsu",    /*   196  gravity wave stress e-w                     */
   "gwd",     /*   197  gravity wave diss                           */
   "src",     /*   198  skin resevoir content                       */
   "sfvc",    /*   199  sfc vegetation cover                        */
   "orgvsg",  /*   200  orgographic variance subgrid                */
   "tasmx",   /*   201  max sfc temp                                */
   "tasmn",   /*   202  min sfc temp                                */
   "v203",    /*   203                                              */
   "praw",    /*   204  precip analysis weights                     */
   "mrro",    /*   205  runoff                                      */
   "cvzz",    /*   206  zz variance                                 */
   "cvtz",    /*   207  tz covariance                               */
   "cvtt",    /*   208  tt variance                                 */
   "cvqz",    /*   209  qz covariance                               */
   "cvqt",    /*   210  qt covariance                               */
   "cvqq",    /*   211  qq variance                                 */
   "cvuz",    /*   212  uz covariance                               */
   "cvut",    /*   213  ut covariance                               */
   "cvuq",    /*   214  uq covariance                               */
   "cvuu",    /*   215  uu variance                                 */
   "cvvz",    /*   216  vz covariance                               */
   "cvvt",    /*   217  vt covariance                               */
   "cvvq",    /*   218  vq covariance                               */
   "cvvu",    /*   219  vu covariance                               */
   "cvvv",    /*   220  vv variance                                 */
   "cvwz",    /*   221  wz covariance                               */
   "cvwt",    /*   222  wt covariance                               */
   "cvwq",    /*   223  wq covariance                               */
   "cvwu",    /*   224  wu covariance                               */
   "cvwv",    /*   225  wv covariance                               */
   "cvww",    /*   226  ww variance                                 */
   "cvrr",    /*   227  rh variance                                 */
   "pr",      /*   228  total precip                                */
   "tauui",   /*   229  instanteous sfc stress u [Nm-2]             */
   "tauvi",   /*   230  instanteous sfc stress v [Nm-2]             */
   "hfssi",   /*   231  instanteous sfc sensible heat flux [Wm-2]   */
   "hflsi",   /*   232  instanteous sfc latent heat flux            */
   "husa",    /*   233  apparent sfc humidity                       */
   "logsfr",  /*   234  log sfc roughness                           */
   "tgs",     /*   235  skin temperature [K]                        */
   "tso4",    /*   236  soil temperature level 4 [K]                */
   "mrso4",   /*   237  soil wetness level 4 [m (H2O)]              */
   "tgs",     /*   238  t of snow layer [K]                         */
   "prsc",    /*   239  convective snow [m]                         */
   "prsl",    /*   240  large scale snow [m]                        */
   "cllw",    /*   241  cloud liquid water                          */
   "clct",    /*   242  total cloud cover                           */
   "albsf",   /*   243  forecast albedo                             */
   "sfrf",    /*   244  forecast sfc roughness                      */
   "logsfcrf",/*   245  log offorecast sfc roughness                */
   "wspds",   /*   246  10 m wind speed                             */   
   "taum",    /*   247  magnitude of momentum flux                  */
   "v248",    /*   248                                              */
   "gwmf",    /*   249  gravity wave dray momentum flux [Nm-2]      */
   "v250",    /*   250                                              */
   "v251",    /*   251                                              */
   "v252",    /*   252                                              */
   "v253",    /*   253                                              */
   "v254",    /*   254                                              */
   "v255"     /*   255                                              */
};



/*
 * EC_ext	v1.0 wesley ebisuzaki
 *
 * prints something readable from the EC stream parameter\
 *
 * prefix and suffix are only printed if EC_ext has text
 */

void EC_ext(unsigned char *pds, char *prefix, char *suffix) {

    /* int i;
    printf("\n");
    for (i=0; i < PDS_LEN(pds); i++) {
      printf("%x ",pds[i]);
    }
    */
    if (PDS_Center(pds) == ECMWF && PDS_LEN(pds) >= 45) {

        switch(PDS_EcStream(pds)) {
	    case 1043:
		printf("%smon mean%s", prefix, suffix);
		break;
	    case 1070:
		printf("%smon (co)var%s", prefix, suffix);
		break;
	    case 1071:
		printf("%smon mean from daily%s", prefix, suffix);
		break;
	    default:
		printf("%sECMWF stream?%s", prefix, suffix);
		break;
	}
    }
}
