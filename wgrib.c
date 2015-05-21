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

/* search order for parameter names
 *
 * #define P_TABLE_FIRST
 * look at external parameter table first
 *
 * otherwise use builtin NCEP-2 or ECMWF-160 first
 */
/* #define P_TABLE_FIRST */

/* search order for external parameter table
 * 1) environment variable GRIBTAB
 * 2) environment variable gribtab
 * 3) the file 'gribtab' in current directory
 */

/* cnames.c */
/* then default values */
char *k5toa(unsigned char *pds);
char *k5_comments(unsigned char *pds);
int setup_user_table(int center, int subcenter, int ptable);


struct ParmTable {
	char *name, *comment;
};

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

int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny);
void GDS_prt_thin_lon(unsigned char *gds);
/* version 3.4 of grib headers  w. ebisuzaki */
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

/* this requires a 32-bit default integer machine */
#define PDS_Field(pds)		((pds[8]<<24)+(pds[9]<<16)+(pds[10]<<8)+pds[11])

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
#define PDS_Year4(pds)          (pds[12] + 100*(pds[24] - (pds[12] != 0)))

/* various centers */
#define NMC			7
#define ECMWF			98

/* ECMWF Extensions */

#define PDS_EcLocalId(pds)	(PDS_LEN(pds) >= 41 ? (pds[40]) : 0)
#define PDS_EcClass(pds)	(PDS_LEN(pds) >= 42 ? (pds[41]) : 0)
#define PDS_EcType(pds)		(PDS_LEN(pds) >= 43 ? (pds[42]) : 0)
#define PDS_EcStream(pds)	(PDS_LEN(pds) >= 45 ? (INT2(pds[43], pds[44])) : 0)



#define VERSION "wgrib v1.5.0b3 (11-21-96) Wesley Ebisuzaki"

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
 * 6/96  - v1.4.1a: faster grib->binary decode, updated ncep parameter table, mod. in clim. desc
 * 7/96  - v1.5.0: parameter-table aware, -v option changed, added "comments"
 *                 increased NTRY to 100 in seek_grib
 * 11/96 - v1.5.0b: added ECMWF parameter table 128
 * 1/97 - v1.5.0b2: if nxny != nx*ny { nx = nxny; ny = 1 }
 *
 */

/*
 * MSEEK = I/O buffer size for seek_grib
 */

#define MSEEK 1024
#define BUFF_ALLOC0	40000


#ifndef min
#define min(a,b)  ((a) < (b) ? (a) : (b))
#define max(a,b)  ((a) < (b) ? (b) : (a))
#endif

int main(int argc, char **argv) {

    unsigned char *buffer;
    float *array;
    double temp, rmin, rmax;
    int i, nx, ny, file_arg;
    long int len_grib, pos = 0, nxny, buffer_size, n_dump, count = 1;
    unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;
    FILE *input, *dump_file = NULL;
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
	fprintf(stderr, "  -v  verbose (verbose inventory)\n");
	fprintf(stderr, "  -V  even more verbose (not a inventory)\n");
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
	fprintf(stderr, "     will scan without -d, -p, or -i options\n");
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
	    GDS_grid(gds, &nx, &ny, &nxny);
	}
	else if (bms != NULL) {
	    nxny = nx = BMS_nxny(bms);
	    ny = 1;
	}
	else {
	    if (BDS_NumBits(bds) == 0) {
                nxny = nx = 1;
                fprintf(stderr,"Missing GDS, constant record .. cannot "
                    "determine number of data points\n");
	    }
	    else {
	        nxny = nx = BDS_NValues(bds);
	    }
	    ny = 1;
	}

#ifdef CHECK_GRIB
        if (BDS_NumBits(bds) != 0) {
	    i = BDS_NValues(bds);
	    if (bms != NULL) {
		i += missing_points(BMS_bitmap(bms),nx*ny);
	    }
	    if (i != nxny) {
	        fprintf(stderr,"grib header at record %ld: two values of nxny %ld %d\n",
			count,nxny,i);
	        fprintf(stderr,"   LEN %d DataStart %d UnusedBits %d #Bits %d\n",
		    BDS_LEN(bds), BDS_DataStart(bds),BDS_UnusedBits(bds),
		    BDS_NumBits(bds));
                return_code = 15;
		nxny = nx = i;
		ny = 1;
	    }
	}
#endif
 
        if (verbose <= 0) {
	    printf("%ld:%ld:d=%2.2d%2.2d%2.2d%2.2d:%s:", count, pos, 
                PDS_Year(pds), PDS_Month(pds), PDS_Day(pds), PDS_Hour(pds),
	        k5toa(pds));
            if (verbose == 0) printf("kpds5=%d:kpds6=%d:kpds7=%d:TR=%d:P1=%d:P2=%d:TimeU=%d:",
	        PDS_PARAM(pds),PDS_KPDS6(pds),PDS_KPDS7(pds),
	        PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                PDS_ForecastTimeUnit(pds));
	    levels(PDS_KPDS6(pds), PDS_KPDS7(pds)); printf(":");
	    PDStimes(PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                PDS_ForecastTimeUnit(pds));
	    if (PDS_Center(pds) == ECMWF) EC_ext(pds,"",":");
	    printf("NAve=%d\n",PDS_NumAve(pds));
       }
       else if (verbose == 1) {
	    printf("%ld:%ld:D=%4.4d%2.2d%2.2d%2.2d:%s:", count, pos, 
                PDS_Year4(pds), PDS_Month(pds), PDS_Day(pds), PDS_Hour(pds),
	        k5toa(pds));
	    levels(PDS_KPDS6(pds), PDS_KPDS7(pds)); printf(":");
            printf("kpds=%d,%d,%d:",
	        PDS_PARAM(pds),PDS_KPDS6(pds),PDS_KPDS7(pds));
	    PDStimes(PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                PDS_ForecastTimeUnit(pds));
	    if (PDS_Center(pds) == ECMWF) EC_ext(pds,"",":");
/*
            printf("grid %d nx %d ny %d :", PDS_Grid(pds),nx,ny);
 */
            printf("\"%s\n", k5_comments(pds));
	}
        else if (verbose == 2) {
	    printf("rec %ld:%ld:date %4.4d%2.2d%2.2d%2.2d %s kpds5=%d "
		"kpds6=%d kpds7=%d levels=(%d,%d) grid=%d ", count, pos, 
                PDS_Year4(pds), PDS_Month(pds), PDS_Day(pds), PDS_Hour(pds),
	        k5toa(pds),
                PDS_PARAM(pds), PDS_KPDS6(pds), PDS_KPDS7(pds), 
                PDS_LEVEL1(pds), PDS_LEVEL2(pds), PDS_Grid(pds));
	    levels(PDS_KPDS6(pds),PDS_KPDS7(pds));
	    printf(" ");
	    if (PDS_Center(pds) == ECMWF) EC_ext(pds,""," ");
	    PDStimes(PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds),
                 PDS_ForecastTimeUnit(pds));
	    if (bms != NULL) 
		printf(" bitmap missing %d", missing_points(BMS_bitmap(bms),nx*ny));
	    printf("\n");
            printf("  %s=%s\n", k5toa(pds), k5_comments(pds));
	
            printf("  timerange %d P1 %d P2 %d TimeU %d  nx %d ny %d GDS grid %d "
		"num_in_ave %d missing %d\n", 
	        PDS_TimeRange(pds),PDS_P1(pds),PDS_P2(pds), 
                PDS_ForecastTimeUnit(pds), nx, ny, 
                gds == NULL ? -1 : GDS_DataType(gds), 
                PDS_NumAve(pds), PDS_NumMissing(pds));

	    printf("  center %d subcenter %d process %d Table %d\n", 
		PDS_Center(pds),PDS_Subcenter(pds),PDS_Model(pds),
                PDS_Vsn(pds));
	    if (gds && GDS_LatLon(gds)) 
		printf("  latlon: lat  %f to %f by %f  nxny %ld\n"
                       "          long %f to %f by %f, (%d x %d) scan %d"
			" bdsgrid %d\n",
		  0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		  0.001*GDS_LatLon_dy(gds), nxny, 0.001*GDS_LatLon_Lo1(gds),
		  0.001*GDS_LatLon_Lo2(gds), 0.001*GDS_LatLon_dx(gds),
	    	  nx, ny, GDS_LatLon_Scan(gds),
		  BDS_Grid(bds));
	    if (gds && GDS_Gaussian(gds) && nx != -1)
		printf("  gaussian: lat  %f to %f\n"
                       "            long %f to %f by %f, (%d x %d) scan %d"
			" bdsgrid %d\n",
		  0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		  0.001*GDS_LatLon_Lo1(gds), 0.001*GDS_LatLon_Lo2(gds), 
		  0.001*GDS_LatLon_dx(gds),
	    	  nx, ny, GDS_LatLon_Scan(gds),
		  BDS_Grid(bds));
	    if (gds && GDS_Gaussian(gds) && nx == -1) {
		printf("  thinned gaussian: lat  %f to %f\n"
                       "     lon %f   %ld grid pts   (%d x %d) scan %d"
			" bdsgrid %d  nlat:\n",
		  0.001*GDS_LatLon_La1(gds), 0.001*GDS_LatLon_La2(gds),
		  0.001*GDS_LatLon_Lo1(gds),
	    	  nxny, nx, ny, GDS_LatLon_Scan(gds),
		  BDS_Grid(bds));
		  GDS_prt_thin_lon(gds);
	    }
	    if (gds && GDS_Polar(gds))
		printf("  polar stereo: Lat1 %f Long1 %f Orient %f\n"
			"     %s pole (%d x %d) Dx %d Dy %d scan %d\n",
		    0.001*GDS_Polar_La1(gds),0.001*GDS_Polar_Lo1(gds),
		    0.001*GDS_Polar_Lov(gds),
		    GDS_Polar_pole(gds) == 0 ? "north" : "south", nx,ny,
		    GDS_Polar_Dx(gds),GDS_Polar_Dy(gds),
		    GDS_Polar_Scan(gds));
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
 * v1.3 8/96 Wesley Ebisuzaki increase NTRY from 3 to 100 for the folks
 *      at Automation decided a 21 byte WMO bulletin header wasn't long 
 *      enough and decided to go to an 8K header.  
 *
 */

#ifndef min
   #define min(a,b)  ((a) < (b) ? (a) : (b))
#endif

#define NTRY 100
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

extern  struct ParmTable parm_table_ncep[256];
extern  struct ParmTable parm_table_omb[256];
extern  struct ParmTable parm_table_ecmwf_128[256];
extern  struct ParmTable parm_table_ecmwf_160[256];
extern  struct ParmTable parm_table_user[256];

char *k5toa(unsigned char *pds) {

    int i, parm, center, subcenter, ptable;
    static int missing_count = 0;

    parm = PDS_PARAM(pds);
    center = PDS_Center(pds);
    subcenter = PDS_Subcenter(pds);
    ptable = PDS_Vsn(pds);

#ifdef P_TABLE_FIRST
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return parm_table_user[parm].name;
    if (center == NMC && ptable <= 3) return parm_table_ncep[parm].name;
    if (center == NMC && ptable == 128) return parm_table_omb[parm].name;
    if (center == ECMWF && ptable == 128) return parm_table_ecmwf_128[parm].name;
    if (center == ECMWF && ptable == 160) return parm_table_ecmwf_160[parm].name;
#else
    if (center == NMC && ptable <= 3) return parm_table_ncep[parm].name;
    if (center == NMC && ptable == 128) return parm_table_omb[parm].name;
    if (center == ECMWF && ptable == 128) return parm_table_ecmwf_128[parm].name;
    if (center == ECMWF && ptable == 160) return parm_table_ecmwf_160[parm].name;
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return parm_table_user[parm].name;
#endif

    if ((ptable > 3 || parm > 127) && missing_count++ == 0) {
	fprintf(stderr,
            "\nUndefined parameter table (center %d-%d table %d), using NCEP-2\n",
            center, subcenter, ptable);
    }
    return parm_table_ncep[parm].name;
}

char *k5_comments(unsigned char *pds) {

    int i, parm, center, subcenter, ptable;

    parm = PDS_PARAM(pds);
    center = PDS_Center(pds);
    subcenter = PDS_Subcenter(pds);
    ptable = PDS_Vsn(pds);

#ifdef P_TABLE_FIRST
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return parm_table_user[parm].comment;
    if (center == NMC && ptable <= 3) return parm_table_ncep[parm].comment;
    if (center == NMC && ptable == 128) return parm_table_omb[parm].comment;
    if (center == ECMWF && ptable == 128) return parm_table_ecmwf_128[parm].comment;
    if (center == ECMWF && ptable == 160) return parm_table_ecmwf_160[parm].comment;
#else
    if (center == NMC && ptable <= 3) return parm_table_ncep[parm].comment;
    if (center == NMC && ptable == 128) return parm_table_omb[parm].comment;
    if (center == ECMWF && ptable == 128) return parm_table_ecmwf_128[parm].comment;
    if (center == ECMWF && ptable == 160) return parm_table_ecmwf_160[parm].comment;
    i = setup_user_table(center, subcenter, ptable);
    if (i == 1) return parm_table_user[parm].comment;
#endif
    return parm_table_ncep[parm].comment;
}


/* 1996				wesley ebisuzaki
 *
 * 4/96 v1.1 faster
 *
 * returns a filled array
 */

static unsigned int mask[] = {0,1,3,7,15,31,63,127,255};

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
	int n_bits, int n, double ref, double scale) {

    int t_bits, c_bits, j_bits;
    unsigned int i, j, map_mask, tbits, jmask, bbits;
    long int jj;

    tbits = bbits = map_mask = 0;

    /* assume integer has 32+ bits */
    if (n_bits <= 25) {
        jmask = (1 << n_bits) - 1;
        t_bits = 0;

        if (bitmap) {
            while (n-- > 0) {
	        if (bitmap) {
		    if (map_mask == 0) {
			map_mask = 128;
			bbits = *bitmap++;
		    }
	            if ((bbits & map_mask) == 0) {
		        *flt++ = UNDEFINED;
	                map_mask >>= 1;
		        continue;
	            }
	            map_mask >>= 1;
	        }
	        while (t_bits < n_bits) {
	            tbits = (tbits * 256) + *bits++;
	            t_bits += 8;
	        }
	        t_bits -= n_bits;
	        j = (tbits >> t_bits) & jmask;
	        *flt++ = ref + scale*j;
            }
        }
        else {
	    for (i = 0; i < n; i++) {
                while (t_bits < n_bits) {
                    tbits = (tbits * 256) + *bits++;
                    t_bits += 8;
                }
                t_bits -= n_bits;
	        j = (tbits >> t_bits) & jmask;
                flt[i] = (tbits >> t_bits) & jmask;
            }
	    /* at least this vectorizes :) */
	    for (i = 0; i < n; i++) {
		flt[i] = ref + scale*flt[i];
	    }
        }
    }
    else {
        c_bits = 8;
        while (n-- > 0) {
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

	    jj = 0;
	    j_bits = n_bits;
	    while (c_bits <= j_bits) {
	        if (c_bits == 8) {
		    jj = (jj << 8) + *bits++;
		    j_bits -= 8;
	        }
	        else {
		    jj = (jj << c_bits) + (*bits & mask[c_bits]);
		    bits++;
		    j_bits -= c_bits;
		    c_bits = 8;
	        }
	    }
	    if (j_bits) {
	        c_bits -= j_bits;
	        jj = (jj << j_bits) + ((*bits >> c_bits) & mask[j_bits]);
	    }
	    *flt++ = ref + scale*jj;
        }
    }
    return;
}

/*
 * convert a float to an ieee single precision number v1.1
 * (big endian)
 *                      Wesley Ebisuzaki
 *
 * bugs: doesn't handle subnormal numbers
 * bugs: assumes length of integer >= 25 bits
 */

int flt2ieee(float x, unsigned char *ieee) {

	int sign, exp;
        unsigned int umant;
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

        /* 2^24 = 16777216 */

	umant = mant * 16777216 + 0.5;
	if (umant >= 16777216) {
            umant = umant / 2;
            exp++;
        }
        /* bit 24 should be a 1 .. not used in ieee format */

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

        ieee[3] = umant & 255;
        ieee[2] = (umant >> 8) & 255;
        ieee[1] = ((exp & 1) << 7) + ((umant >> 16) & 127);
	return 0;
}


/* wesley ebisuzaki v1.0
 *
 * write ieee file -- big endian format
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
	int fcst_len = 0;

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
		if (p1 != 0) printf("start@%d%s:",p1,unit);
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
	case 51: if (p1 == 0) {
		    printf("clim %d%s:",p2,unit);
		}
		else if (p1 == 1) {
		    printf("clim (diurnal) %d%s:",p2,unit);
		}
		else {
		    printf("clim? p1=%d? %d%s?:",p1,p2,unit);
		}
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
 *  number of missing data points w. ebisuzaki
 *
 *  v1.1: just faster my dear
 *  v1.2: just faster my dear
 *
 */

static int bitsum[256] = {
    8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0};


int missing_points(unsigned char *bitmap, int n) {

    int count;
    unsigned int tmp;
    if (bitmap == NULL) return 0;

    count = 0;
    while (n >= 8) {
	tmp = *bitmap++;
	n -= 8;
        count += bitsum[tmp];
    }
    tmp = *bitmap | ((1 << (8 - n)) - 1);
    count += bitsum[tmp];

    return count;
}

struct ParmTable parm_table_ncep[256] = {
  {"var0", "undefined"}, /* 0 */
  {"PRES", "Pressure [Pa]"},     /* 1 */
  {"PRMSL", "Pressure reduced to MSL [Pa]"},     /* 2 */
  {"PTEND", "Pressure tendency [Pa/s]"},     /* 3 */
  {"var4", "undefined"}, /* 4 */
  {"var5", "undefined"}, /* 5 */
  {"GP", "Geopotential [m^2/s^2]"},     /* 6 */
  {"HGT", "Geopotential height [gpm]"},     /* 7 */
  {"DIST", "Geometric height [m]"},     /* 8 */
  {"HSTDV", "Std dev of height [m]"},     /* 9 */
  {"HVAR", "Variance of height [m^2]"},     /* 10 */
  {"TMP", "Temp. [K]"},     /* 11 */
  {"VTMP", "Virtual temp. [K]"},     /* 12 */
  {"POT", "Potential temp. [K]"},     /* 13 */
  {"EPOT", "Pseudo-adiabatic pot. temp. [K]"},     /* 14 */
  {"TMAX", "Max. temp. [K]"},     /* 15 */
  {"TMIN", "Min. temp. [K]"},     /* 16 */
  {"DPT", "Dew point temp. [K]"},     /* 17 */
  {"DEPR", "Dew point depression [K]"},     /* 18 */
  {"LAPR", "Lapse rate [K/m]"},     /* 19 */
  {"VISIB", "Visibility [m]"},     /* 20 */
  {"RDSP1", "Radar spectra (1) [non-dim]"},     /* 21 */
  {"RDSP2", "Radar spectra (2) [non-dim]"},     /* 22 */
  {"RDSP3", "Radar spectra (3) [non-dim]"},     /* 23 */
  {"var24", "undefined"}, /* 24 */
  {"TMPA", "Temp. anomaly [K]"},     /* 25 */
  {"PRESA", "Pressure anomaly [Pa]"},     /* 26 */
  {"GPA", "Geopotential height anomaly [gpm]"},     /* 27 */
  {"WVSP1", "Wave spectra (1) [non-dim]"},     /* 28 */
  {"WVSP2", "Wave spectra (2) [non-dim]"},     /* 29 */
  {"WVSP3", "Wave spectra (3) [non-dim]"},     /* 30 */
  {"WDIR", "Wind direction [degree]"},     /* 31 */
  {"WIND", "Wind speed [m/s]"},     /* 32 */
  {"UGRD", "u wind [m/s]"},     /* 33 */
  {"VGRD", "v wind [m/s]"},     /* 34 */
  {"STRM", "Stream function [m^2/s]"},     /* 35 */
  {"VPOT", "Velocity potential [m^2/s]"},     /* 36 */
  {"MNTSF", "Montgomery stream function [m^2/s^2]"},     /* 37 */
  {"SGCVV", "Sigma coord. vertical velocity [/s]"},     /* 38 */
  {"VVEL", "Pressure vertical velocity [Pa/s]"},     /* 39 */
  {"DZDT", "Geometric vertical velocity [m/s]"},     /* 40 */
  {"ABSV", "Absolute vorticity [/s]"},     /* 41 */
  {"ABSD", "Absolute divergence [/s]"},     /* 42 */
  {"RELV", "Relative vorticity [/s]"},     /* 43 */
  {"RELD", "Relative divergence [/s]"},     /* 44 */
  {"VUCSH", "Vertical u shear [/s]"},     /* 45 */
  {"VVCSH", "Vertical v shear [/s]"},     /* 46 */
  {"DIRC", "Direction of current [degree]"},     /* 47 */
  {"SPC", "Speed of current [m/s]"},     /* 48 */
  {"UOGRD", "u of current [m/s]"},     /* 49 */
  {"VOGRD", "v of current [m/s]"},     /* 50 */
  {"SPFH", "Specific humidity [kg/kg]"},     /* 51 */
  {"RH", "Relative humidity [%]"},     /* 52 */
  {"MIXR", "Humidity mixing ratio [kg/kg]"},     /* 53 */
  {"PWAT", "Precipitable water [kg/m^2]"},     /* 54 */
  {"VAPP", "Vapor pressure [Pa]"},     /* 55 */
  {"SATD", "Saturation deficit [Pa]"},     /* 56 */
  {"EVP", "Evaporation [kg/m^2]"},     /* 57 */
  {"CICE", "Cloud Ice [kg/m^2]"},     /* 58 */
  {"PRATE", "Precipitation rate [kg/m^2/s]"},     /* 59 */
  {"TSTM", "Thunderstorm probability [%]"},     /* 60 */
  {"APCP", "Total precipitation [kg/m^2]"},     /* 61 */
  {"NCPCP", "Large scale precipitation [kg/m^2]"},     /* 62 */
  {"ACPCP", "Convective precipitation [kg/m^2]"},     /* 63 */
  {"SRWEQ", "Snowfall rate water equiv. [kg/m^2/s]"},     /* 64 */
  {"WEASD", "Accum. snow [kg/m^2]"},     /* 65 */
  {"SNOD", "Snow depth [m]"},     /* 66 */
  {"MIXHT", "Mixed layer depth [m]"},     /* 67 */
  {"TTHDP", "Transient thermocline depth [m]"},     /* 68 */
  {"MTHD", "Main thermocline depth [m]"},     /* 69 */
  {"MTHA", "Main thermocline anomaly [m]"},     /* 70 */
  {"TCDC", "Total cloud cover [%]"},     /* 71 */
  {"CDCON", "Convective cloud cover [%]"},     /* 72 */
  {"LCDC", "Low level cloud cover [%]"},     /* 73 */
  {"MCDC", "Mid level cloud cover [%]"},     /* 74 */
  {"HCDC", "High level cloud cover [%]"},     /* 75 */
  {"CWAT", "Cloud water [kg/m^2]"},     /* 76 */
  {"var77", "undefined"}, /* 77 */
  {"SNOC", "Convective snow [kg/m^2]"},     /* 78 */
  {"SNOL", "Large scale snow [kg/m^2]"},     /* 79 */
  {"WTMP", "Water temp. [K]"},     /* 80 */
  {"LAND", "Land-sea mask [1=land; 0=sea]"},     /* 81 */
  {"DSLM", "Deviation of sea level from mean [m]"},     /* 82 */
  {"SFCR", "Surface roughness [m]"},     /* 83 */
  {"ALBDO", "Albedo [%]"},     /* 84 */
  {"TSOIL", "Soil temp. [K]"},     /* 85 */
  {"SOILM", "Soil moisture content [kg/m^2]"},     /* 86 */
  {"VEG", "Vegetation [%]"},     /* 87 */
  {"SALTY", "Salinity [kg/kg]"},     /* 88 */
  {"DEN", "Density [kg/m^2]"},     /* 89 */
  {"RUNOF", "Runoff [kg/m^2]"},     /* 90 */
  {"ICEC", "Ice concentration [ice=1; no ice=0]"},     /* 91 */
  {"ICETK", "Ice thickness [m]"},     /* 92 */
  {"DICED", "Direction of ice drift [degree]"},     /* 93 */
  {"SICED", "Speed of ice drift [m/s]"},     /* 94 */
  {"UICE", "u of ice drift [m/s]"},     /* 95 */
  {"VICE", "v of ice drift [m/s]"},     /* 96 */
  {"ICEG", "Ice growth [m]"},     /* 97 */
  {"ICED", "Ice divergence [/s]"},     /* 98 */
  {"SNOM", "Snow melt [kg/m^2]"},     /* 99 */
  {"HTSGW", "Sig height of wind waves and swell [m]"},     /* 100 */
  {"WVDIR", "Direction of wind waves [degree]"},     /* 101 */
  {"WVHGT", "Sig height of wind waves [m]"},     /* 102 */
  {"WVPER", "Mean period of wind waves [s]"},     /* 103 */
  {"SWDIR", "Direction of swell waves [degree]"},     /* 104 */
  {"SWELL", "Sig height of swell waves [m]"},     /* 105 */
  {"SWPER", "Mean period of swell waves [s]"},     /* 106 */
  {"DIRPW", "Primary wave direction [degree]"},     /* 107 */
  {"PERPW", "Primary wave mean period [s]"},     /* 108 */
  {"DIRSW", "Secondary wave direction [degree]"},     /* 109 */
  {"PERSW", "Secondary wave mean period [s]"},     /* 110 */
  {"NSWRS", "Net short wave (surface) [W/m^2]"},     /* 111 */
  {"NLWRS", "Net long wave (surface) [W/m^2]"},     /* 112 */
  {"NSWRT", "Net short wave (top) [W/m^2]"},     /* 113 */
  {"NLWRT", "Net long wave (top) [W/m^2]"},     /* 114 */
  {"LWAVR", "Long wave [W/m^2]"},     /* 115 */
  {"SWAVR", "Short wave [W/m^2]"},     /* 116 */
  {"GRAD", "Global radiation [W/m^2]"},     /* 117 */
  {"var118", "undefined"}, /* 118 */
  {"var119", "undefined"}, /* 119 */
  {"var120", "undefined"}, /* 120 */
  {"LHTFL", "Latent heat flux [W/m^2]"},     /* 121 */
  {"SHTFL", "Sensible heat flux [W/m^2]"},     /* 122 */
  {"BLYDP", "Boundary layer dissipation [W/m^2]"},     /* 123 */
  {"UFLX", "Zonal momentum flux [N/m^2]"},     /* 124 */
  {"VFLX", "Meridional momentum flux [N/m^2]"},     /* 125 */
  {"WMIXE", "Wind mixing energy [J]"},     /* 126 */
  {"IMGD", "Image data [integer]"},     /* 127 */
  {"MSLSA", "Mean sea level pressure (Std Atm) [Pa]"},     /* 128 */
  {"MSLMA", "Mean sea level pressure (MAPS) [Pa]"},     /* 129 */
  {"MSLET", "Mean sea level pressure (ETA model) [Pa]"},     /* 130 */
  {"LFTX", "Surface lifted index [K]"},     /* 131 */
  {"4LFTX", "Best (4-layer) lifted index [K]"},     /* 132 */
  {"KX", "K index [K]"},     /* 133 */
  {"SX", "Sweat index [K]"},     /* 134 */
  {"MCONV", "Horizontal moisture divergence [kg/kg/s]"},     /* 135 */
  {"VSSH", "Vertical speed shear [1/s]"},     /* 136 */
  {"TSLSA", "3-hr pressure tendency [Pa/s]"},     /* 137 */
  {"BVF2", "Brunt-Vaisala frequency^2 [1/s^2]"},     /* 138 */
  {"PVMW", "Potential vorticity (mass-weighted) [1/s/m]"},     /* 139 */
  {"CRAIN", "Categorical rain [yes=1;no=0]"},     /* 140 */
  {"CFRZR", "Categorical freezing rain [yes=1;no=0]"},     /* 141 */
  {"CICEP", "Categorical ice pellets [yes=1;no=0]"},     /* 142 */
  {"CSNOW", "Categorical snow [yes=1;no=0]"},     /* 143 */
  {"SOILW", "Volumetric soil moisture [fraction]"},     /* 144 */
  {"PEVPR", "Potential evaporation rate [w/m^2]"},     /* 145 */
  {"CWORK", "Cloud work function [J/Kg]"},     /* 146 */
  {"U-GWD", "Zonal gravity wave stress [N/m^2]"},     /* 147 */
  {"V-GWD", "Meridional gravity wave stress [N/m^2]"},     /* 148 */
  {"PV___", "Potential vorticity [m^2/s/kg]"},     /* 149 */
  {"var150", "undefined"}, /* 150 */
  {"var151", "undefined"}, /* 151 */
  {"var152", "undefined"}, /* 152 */
  {"MFXDV", "Moisture flux divergence [gr/gr*m/s/m]"},     /* 153 */
  {"var154", "undefined"}, /* 154 */
  {"GFLUX", "Ground heat flux [W/m^2]"},     /* 155 */
  {"CIN", "Convective inhibition [J/kg]"},     /* 156 */
  {"CAPE", "Convective Avail. Pot. Energy [J/kg]"},     /* 157 */
  {"TKE", "Turbulent kinetic energy [J/kg]"},     /* 158 */
  {"CONDP", "Lifted parcel condensation pressure [Pa]"},     /* 159 */
  {"CSUSF", "Clear sky upward solar flux [W/m^2]"},     /* 160 */
  {"CSDSF", "Clear sky downward solar flux [W/m^2]"},     /* 161 */
  {"CSULF", "Clear sky upward long wave flux [W/m^2]"},     /* 162 */
  {"CSDLF", "Clear sky downward long wave flux [W/m^2]"},     /* 163 */
  {"CFNSF", "Cloud forcing net solar flux [W/m^2]"},     /* 164 */
  {"CFNLF", "Cloud forcing net long wave flux [W/m^2]"},     /* 165 */
  {"VBDSF", "Visible beam downward solar flux [W/m^2]"},     /* 166 */
  {"VDDSF", "Visible diffuse downward solar flux [W/m^2]"},     /* 167 */
  {"NBDSF", "Near IR beam downward solar flux [W/m^2]"},     /* 168 */
  {"NDDSF", "Near IR diffuse downward solar flux [W/m^2]"},     /* 169 */
  {"USTR", "U wind stress [N/m^2]"},     /* 170 */
  {"VSTR", "V wind stress [N/m^2]"},     /* 171 */
  {"MFLX", "Momentum flux [N/m^2]"},     /* 172 */
  {"LMH", "Mass point model surface [integer]"},     /* 173 */
  {"LMV", "Velocity point model surface [integer]"},     /* 174 */
  {"SGLYR", "Nearby model level [integer]"},     /* 175 */
  {"NLAT", "Latitude [degree]"},     /* 176 */
  {"NLON", "Longitude [degree]"},     /* 177 */
  {"UMAS", "Mass weighted u [gm/m*K*s]"},     /* 178 */
  {"VMAS", "Mass weighted v [gm/m*K*s]"},     /* 179 */
  {"XPRATE", "corrected precip [Kg/m^2/s]"}, /* 180 */
  {"LPSX", "x-gradient of log pressure [1/m]"},     /* 181 */
  {"LPSY", "y-gradient of log pressure [1/m]"},     /* 182 */
  {"HGTX", "x-gradient of height [m/m]"},     /* 183 */
  {"HGTY", "y-gradient of height [m/m]"},     /* 184 */
  {"STDZ", "Std dev of Geop. hgt. [m]"},     /* 185 */
  {"STDU", "Std dev of zonal wind [m/s]"},     /* 186 */
  {"STDV", "Std dev of meridional wind [m/s]"},     /* 187 */
  {"STDQ", "Std dev of spec. hum. [gm/gm]"},     /* 188 */
  {"STDT", "Std dev of temp. [K]"},     /* 189 */
  {"CBUW", "Covar. u and omega [m/s*Pa/s]"},     /* 190 */
  {"CBVW", "Covar. v and omega [m/s*Pa/s]"},     /* 191 */
  {"CBUQ", "Covar. u and specific hum [m/s*gm/gm]"},     /* 192 */
  {"CBVQ", "Covar. v and specific hum [m/s*gm/gm]"},     /* 193 */
  {"CBTW", "Covar. T and omega [K*Pa/s]"},     /* 194 */
  {"CBQW", "Covar. spec. hum and omega [gm/gm*Pa/s]"},     /* 195 */
  {"CBMZW", "Covar. v and u [m^2/s^2]"},     /* 196 */
  {"CBTZW", "Covar. u and T [K*m/s]"},     /* 197 */
  {"CBTMW", "Covar. v and T [K*m/s]"},     /* 198 */
  {"STDRH", "Std dev of Rel. Hum. [%]"},     /* 199 */
  {"SDTZ", "Std dev of time tend of geop. hgt [m]"},     /* 200 */
  {"ICWAT", "Ice-free water surface [%]"},     /* 201 */
  {"SDTU", "Std dev of time tend of zonal wind [m/s]"},     /* 202 */
  {"SDTV", "Std dev of time tend of merid wind [m/s]"},     /* 203 */
  {"DSWRF", "Downward solar radiation flux [W/m^2]"},     /* 204 */
  {"DLWRF", "Downward long wave flux [W/m^2]"},     /* 205 */
  {"SDTQ", "Std dev of time tend of spec. hum [gm/gm]"},     /* 206 */
  {"MSTAV", "Moisture availability [%]"},     /* 207 */
  {"SFEXC", "Exchange coefficient [kg*m/m^3/s]"},     /* 208 */
  {"MIXLY", "Number of mixed layers next to sfc [integer]"},     /* 209 */
  {"SDTT", "Std dev of time tend of temp. [K]"},     /* 210 */
  {"USWRF", "Upward solar radiation flux [W/m^2]"},     /* 211 */
  {"ULWRF", "Upward long wave flux [W/m^2]"},     /* 212 */
  {"CDLYR", "Non-convective cloud [%]"},     /* 213 */
  {"CPRAT", "Convective precip. rate [kg/m^2/s]"},     /* 214 */
  {"TTDIA", "Temp. tendency by all physics [K/s]"},     /* 215 */
  {"TTRAD", "Temp. tendency by all radiation [K/s]"},     /* 216 */
  {"TTPHY", "Temp. tendency by nonrad physics [K/s]"},     /* 217 */
  {"PREIX", "Precipitation index [fraction]"},     /* 218 */
  {"TSD1D", "Std dev of IR T over 1x1 deg area [K]"},     /* 219 */
  {"NLSGP", "Natural log of surface pressure [ln(kPa)]"},     /* 220 */
  {"SDTRH", "Std dev of time tend of rel hum [%]"},     /* 221 */
  {"5WAVH", "5-wave geopotential height [gpm]"},     /* 222 */
  {"CWAT", "Plant canopy surface water [kg/m^2]"},     /* 223 */
  {"PLTRS", "Max. stomato plant resistance [s/m]"},     /* 224 */
  {"RHCLD", "RH-type cloud cover [%]"},     /* 225 */
  {"BMIXL", "Blackadar's mixing length scale [m]"},     /* 226 */
  {"AMIXL", "Asymptotic mixing length scale [m]"},     /* 227 */
  {"PEVAP", "Pot. evaporation [kg/m^2]"},     /* 228 */
  {"SNOHF", "Snow melt heat flux [W/m^2]"},     /* 229 */
  {"SNOEV", "Snow sublimation heat flux [W/m^2]"},     /* 230 */
  {"MFLUX", "Convective cloud mass flux [Pa/s]"},     /* 231 */
  {"DTRF", "Downward total radiation flux [W/m^2]"},     /* 232 */
  {"UTRF", "Upward total radiation flux [W/m^2]"},     /* 233 */
  {"BGRUN", "Baseflow-groundwater runoff [kg/m^2]"},     /* 234 */
  {"SSRUN", "Storm surface runoff [kg/m^2]"},     /* 235 */
  {"var236", "undefined"}, /* 236 */
  {"OZONE", "Total column ozone [Dobson]"},     /* 237 */
  {"SNOC", "Snow cover [%]"},     /* 238 */
  {"SNOT", "Snow temp. [K]"},     /* 239 */
  {"GLCR", "Permanent snow points [mask]"},     /* 240 */
  {"LRGHR", "Large scale condensation heating [K/s]"},     /* 241 */
  {"CNVHR", "Deep convective heating [K/s]"},     /* 242 */
  {"CNVMR", "Deep convective moistening rate [kg/kg/s]"},     /* 243 */
  {"SHAHR", "Shallow convective heating [K/s]"},     /* 244 */
  {"SHAMR", "Shallow convective moistening rate [kg/kg/s]"},     /* 245 */
  {"VDFHR", "Vertical diffusion heating [K/s]"},     /* 246 */
  {"VDFUA", "Vertical diffusion zonal accel [m/s/s]"},     /* 247 */
  {"VDFVA", "Vertical diffusion meridional accel [m/s/s]"},     /* 248 */
  {"VDFMR", "Vertical diffusion moistening rate [kg/kg/s]"},     /* 249 */
  {"SWHR", "Solar radiative heating [K/s]"},     /* 250 */
  {"LWHR", "Longwave radiative heating [K/s]"},     /* 251 */
  {"CD", "Drag coefficient [non-dim]"},     /* 252 */
  {"FRICV", "Friction velocity [m/s]"},     /* 253 */
  {"RI", "Richardson number [non-dim]"},     /* 254 */
  {"var255", "undefined"}, /* 255 */
};

struct ParmTable parm_table_ecmwf_128[256] = {
  {"var0", "undefined"}, /* 0 */
  {"PRES", "Pressure [Pa]"},     /* 1 */
  {"PRMSL", "Pressure reduced to MSL [Pa]"},     /* 2 */
  {"PTEND", "Pressure tendency [Pa/s]"},     /* 3 */
  {"var4", "undefined"}, /* 4 */
  {"var5", "undefined"}, /* 5 */
  {"GP", "Geopotential [m**2/s**2]"},     /* 6 */
  {"HGT", "Geopotential height [gpm]"},     /* 7 */
  {"DIST", "Geometric height [m]"},     /* 8 */
  {"HSTDV", "Standard deviation of height [m]"},     /* 9 */
  {"HVAR", "Variance of height [m**2]"},     /* 10 */
  {"TMP", "Temperature [K]"},     /* 11 */
  {"VTMP", "Virtual temperature [K]"},     /* 12 */
  {"POT", "Potential temperature [K]"},     /* 13 */
  {"EPOT", "Pseudo-adiabatic potential temperature [K]"},     /* 14 */
  {"TMAX", "Maximum temperature [K]"},     /* 15 */
  {"TMIN", "Minimum temperature [K]"},     /* 16 */
  {"DPT", "Dew point temperature [K]"},     /* 17 */
  {"DEPR", "Dew point depression [K]"},     /* 18 */
  {"LAPR", "Lapse rate [K/m]"},     /* 19 */
  {"VISIB", "Visibility [m]"},     /* 20 */
  {"RDSP1", "Radar spectra (1) [dimensionless]"},     /* 21 */
  {"RDSP2", "Radar spectra (2) [dimensionless]"},     /* 22 */
  {"RDSP3", "Radar spectra (3) [dimensionless]"},     /* 23 */
  {"var24", "undefined"}, /* 24 */
  {"TMPA", "Temperature anomaly [K]"},     /* 25 */
  {"PRESA", "Pressure anomaly [Pa]"},     /* 26 */
  {"GPA", "Geopotential height anomaly [gpm]"},     /* 27 */
  {"WVSP1", "Wave spectra (1) [dimensionless]"},     /* 28 */
  {"WVSP2", "Wave spectra (2) [dimensionless]"},     /* 29 */
  {"WVSP3", "Wave spectra (3) [dimensionless]"},     /* 30 */
  {"WDIR", "Wind direction [degree]"},     /* 31 */
  {"WIND", "Wind speed [m/s]"},     /* 32 */
  {"UGRD", "u wind [m/s]"},     /* 33 */
  {"VGRD", "v wind [m/s]"},     /* 34 */
  {"STRM", "Stream function [m**2/s]"},     /* 35 */
  {"VPOT", "Velocity potential [m**2/s]"},     /* 36 */
  {"MNTSF", "Montgomery stream function [m**2/s**2]"},     /* 37 */
  {"SGCVV", "Sigma coord. vertical velocity [/s]"},     /* 38 */
  {"VVEL", "Pressure vertical velocity [Pa/s]"},     /* 39 */
  {"DZDT", "Geometric vertical velocity [m/s]"},     /* 40 */
  {"ABSV", "Absolute vorticity [/s]"},     /* 41 */
  {"ABSD", "Absolute divergence [/s]"},     /* 42 */
  {"RELV", "Relative vorticity [/s]"},     /* 43 */
  {"RELD", "Relative divergence [/s]"},     /* 44 */
  {"VUCSH", "Vertical u shear [/s]"},     /* 45 */
  {"VVCSH", "Vertical v shear [/s]"},     /* 46 */
  {"DIRC", "Direction of current [degree]"},     /* 47 */
  {"SPC", "Speed of current [m/s]"},     /* 48 */
  {"UOGRD", "u of current [m/s]"},     /* 49 */
  {"VOGRD", "v of current [m/s]"},     /* 50 */
  {"SPFH", "Specific humidity [kg/kg]"},     /* 51 */
  {"RH", "Relative humidity [percent]"},     /* 52 */
  {"MIXR", "Humidity mixing ratio [kg/kg]"},     /* 53 */
  {"PWAT", "Precipitable water [kg/m**2]"},     /* 54 */
  {"VAPP", "Vapor pressure [Pa]"},     /* 55 */
  {"SATD", "Saturation deficit [Pa]"},     /* 56 */
  {"EVP", "Evaporation [kg/m**2]"},     /* 57 */
  {"CICE", "Cloud Ice [kg/m**2]"},     /* 58 */
  {"PRATE", "Precipitation rate [kg/m**2/s]"},     /* 59 */
  {"TSTM", "Thunderstorm probability [percent]"},     /* 60 */
  {"APCP", "Total precipitation [kg/m**2]"},     /* 61 */
  {"NCPCP", "Large scale precipitation [kg/m**2]"},     /* 62 */
  {"ACPCP", "Convective precipitation [kg/m**2]"},     /* 63 */
  {"SRWEQ", "Snowfall rate water equivalent [kg/m**2/s]"},     /* 64 */
  {"WEASD", "Water equiv. of accum. snow depth [kg/m**2]"},     /* 65 */
  {"SNOD", "Snow depth [m]"},     /* 66 */
  {"MIXHT", "Mixed layer depth [m]"},     /* 67 */
  {"TTHDP", "Transient thermocline depth [m]"},     /* 68 */
  {"MTHD", "Main thermocline depth [m]"},     /* 69 */
  {"MTHA", "Main thermocline anomaly [m]"},     /* 70 */
  {"TCDC", "Total cloud cover [percent]"},     /* 71 */
  {"CDCON", "Convective cloud cover [percent]"},     /* 72 */
  {"LCDC", "Low level cloud cover [percent]"},     /* 73 */
  {"MCDC", "Mid level cloud cover [percent]"},     /* 74 */
  {"HCDC", "High level cloud cover [percent]"},     /* 75 */
  {"CWAT", "Cloud water [kg/m**2]"},     /* 76 */
  {"var77", "undefined"}, /* 77 */
  {"SNOC", "Convective snow [kg/m**2]"},     /* 78 */
  {"SNOL", "Large scale snow [kg/m**2]"},     /* 79 */
  {"WTMP", "Water temperature [K]"},     /* 80 */
  {"LAND", "Land-sea mask (1=land; 0=sea) [integer]"},     /* 81 */
  {"DSLM", "Deviation of sea level from mean [m]"},     /* 82 */
  {"SFCR", "Surface roughness [m]"},     /* 83 */
  {"ALBDO", "Albedo [percent]"},     /* 84 */
  {"TSOIL", "Soil temperature [K]"},     /* 85 */
  {"SOILM", "Soil moisture content [kg/m**2]"},     /* 86 */
  {"VEG", "Vegetation [percent]"},     /* 87 */
  {"SALTY", "Salinity [kg/kg]"},     /* 88 */
  {"DEN", "Density [kg/m**2]"},     /* 89 */
  {"RUNOF", "Runoff [kg/m**2]"},     /* 90 */
  {"ICEC", "Ice concentration (ice=1; no ice=0) [1/0]"},     /* 91 */
  {"ICETK", "Ice thickness [m]"},     /* 92 */
  {"DICED", "Direction of ice drift [degree]"},     /* 93 */
  {"SICED", "Speed of ice drift [m/s]"},     /* 94 */
  {"UICE", "u of ice drift [m/s]"},     /* 95 */
  {"VICE", "v of ice drift [m/s]"},     /* 96 */
  {"ICEG", "Ice growth [m]"},     /* 97 */
  {"ICED", "Ice divergence [/s]"},     /* 98 */
  {"SNOM", "Snow melt [kg/m**2]"},     /* 99 */
  {"HTSGW", "Sig height of wind waves and swell [m]"},     /* 100 */
  {"WVDIR", "Direction of wind waves [degree]"},     /* 101 */
  {"WVHGT", "Significant height of wind waves [m]"},     /* 102 */
  {"WVPER", "Mean period of wind waves [s]"},     /* 103 */
  {"SWDIR", "Direction of swell waves [degree]"},     /* 104 */
  {"SWELL", "Significant height of swell waves [m]"},     /* 105 */
  {"SWPER", "Mean period of swell waves [s]"},     /* 106 */
  {"DIRPW", "Primary wave direction [degree]"},     /* 107 */
  {"PERPW", "Primary wave mean period [s]"},     /* 108 */
  {"DIRSW", "Secondary wave direction [degree]"},     /* 109 */
  {"PERSW", "Secondary wave mean period [s]"},     /* 110 */
  {"NSWRS", "Net short wave radiation (surface) [W/m**2]"},     /* 111 */
  {"NLWRS", "Net long wave radiation (surface) [W/m**2]"},     /* 112 */
  {"NSWRT", "Net short wave radiation (top) [W/m**2]"},     /* 113 */
  {"NLWRT", "Net long wave radiation (top) [W/m**2]"},     /* 114 */
  {"LWAVR", "Long wave radiation [W/m**2]"},     /* 115 */
  {"SWAVR", "Short wave radiation [W/m**2]"},     /* 116 */
  {"GRAD", "Global radiation [W/m**2]"},     /* 117 */
  {"var118", "undefined"}, /* 118 */
  {"var119", "undefined"}, /* 119 */
  {"var120", "undefined"}, /* 120 */
  {"LHTFL", "Latent heat flux [W/m**2]"},     /* 121 */
  {"SHTFL", "Sensible heat flux [W/m**2]"},     /* 122 */
  {"BLYDP", "Boundary layer dissipation [W/m**2]"},     /* 123 */
  {"UFLX", "Zonal component of momentum flux [N/m**2]"},     /* 124 */
  {"VFLX", "Meridional component of momentum flux [N/m**2]"},     /* 125 */
  {"WMIXE", "Wind mixing energy [J]"},     /* 126 */

  /* ECMWF local table */

  {"AT",      /* 127 */ "Atmospheric Tide"},
  {"BV",      /* 128 */ "Budget Values"},
  {"Z",       /* 129 */ "Geopotential [m2 s-2]"},
  {"T",       /* 130 */ "Temperature [K]"},
  {"U",       /* 131 */ "U-component of Wind [ms-1]"},
  {"V",       /* 132 */ "V-component of Wind [ms-1]"},
  {"Q",       /* 133 */ "Specific Humidity [kg/kg]"},
  {"SP",      /* 134 */ "Surface Pressure [Pa]"},
  {"W",       /* 135 */ "Vertical Velocity [Pa s-1]"},
  {"TCW",     /* 136 */ "Total column water (vapor+drops+ice) [kg/m2]"},
  {"TCWV",    /* 137 */ "Total column water vapor [kg/m2]"},
  {"VO",      /* 138 */ "relative vorticity [s-1]"},
  {"STL1",    /* 139 */ "soil temperature level 1 [K]"},
  {"SWL1",    /* 140 */ "soil moisture level 1 [m (H20)]"},
  {"SD",      /* 141 */ "Snow Depth [m]"},
  {"LSP",     /* 142 */ "Large Scale Precipitation [m]"},
  {"CP",      /* 143 */ "Convective Precipitation [m]"},
  {"SF",      /* 144 */ "Snow Fall [m]"},
  {"BLD",     /* 145 */ "Boundary Layer Dissipation [Wm-2]"},
  {"SSHF",    /* 146 */ "Surface Flux of Sensible Heat [Wm-2]"},
  {"SLHF",    /* 147 */ "Surface Flux of Latent Heat [Wm-2]"},
  {"v148",    /* 148 */ "undefined"},
  {"v149",    /* 149 */ "undefined"},
  {"v150",    /* 150 */ "undefined"},
  {"MSL",     /* 151 */ "Mean Sea Level (MSL) Pressure [Pa]    Pa"},
  {"LNSP",    /* 152 */ "Log Surface Pressure"},
  {"v153",    /* 153 */ "undefined"},
  {"v154",    /* 154 */ "undefined"},
  {"D",       /* 155 */ "Divergence [s-1]"},
  {"GH",      /* 156 */ "Height (Geopotential) [m]"},
  {"R",       /* 157 */ "Relative Humidity [%]"},
  {"TSP",     /* 158 */ "Tendency of Surface Pressure [Pa s-1]"},
  {"v159",    /* 159 */ "undefined"},
  {"SDOR",    /* 160 */ "Standard deviation of orography"},
  {"ISOR",    /* 161 */ "Anisotropy of subgrid scale orography"},
  {"ANOR",    /* 162 */ "Angle of subgrid scale orography"},
  {"SLOR",    /* 163 */ "Slope of subgrid scale orography"},
  {"TCC",     /* 164 */ "cloud cover total [0-1]"},
  {"10U",     /* 165 */ "U-wind at 10 [ms-1]"},
  {"10V",     /* 166 */ "V-wind at 10 [ms-1]"},
  {"2T",      /* 167 */ "Temperature at 2 m [K]"},
  {"2D",      /* 168 */ "Dewpoint at 2 m [K]"},
  {"v169",    /* 169 */ "undefined"},
  {"STL2",    /* 170 */ "soil temperature level 2 [K]"},
  {"SWL2",    /* 171 */ "soil wetness level 2 [m (H20)]"},
  {"LSM",     /* 172 */ "land-sea mask [(0,1)]"},
  {"SR",      /* 173 */ "sfc roughness [m]"},
  {"AL",      /* 174 */ "Albedo [0-1]"},
  {"v175",    /* 175 */ "undefined"},
  {"SSR",     /* 176 */ "Net Shortwave Radiation (surface) [Wm-2]"},
  {"STR",     /* 177 */ "Net Longwave Radiation (surface) [Wm-2]"},
  {"TSR",     /* 178 */ "Net Shortwave Radiation (toa) [Wm-2"},
  {"TTR",     /* 179 */ "Net Longwave Radiation (toa) [Wm-2]"},
  {"EWSS",    /* 180 */ "U-component of Surface Wind Stress [Nm-2]"},
  {"NSSS",    /* 181 */ "V-component of Surface Wind Stress [Nm-2]"},
  {"E",       /* 182 */ "Evaporation [m (H2O)]"},
  {"STL3",    /* 183 */ "soil temp level 3 [m (H2O)]"},
  {"SWL3",    /* 184 */ "soil moisture level 3 [K]"},
  {"CCC",     /* 185 */ "cloud convective [0-1]"},
  {"LCC",     /* 186 */ "cloud low [0-1]"},
  {"MCC",     /* 187 */ "cloud mid [0-1]"},
  {"HCC",     /* 188 */ "cloud high [0-1]"},
  {"v189",    /* 189 */ "undefined"},
  {"EWOV",    /* 190 */ "orographic variance e-w [m2]"},
  {"NSOV",    /* 191 */ "orographic variance n-s [m2]"},
  {"NWOV",    /* 192 */ "orographic variance nw-se [m2]"},
  {"NEOV",    /* 193 */ "orographic variance ne-sw [m2]"},
  {"v194",    /* 194 */ "undefined"},
  {"LGWS",    /* 195 */ "gravity wave stress n-s [n/m2s]"},
  {"MGWS",    /* 196 */ "gravity wave stress e-w [n/m2s]"},
  {"GWD",     /* 197 */ "gravity wave diss [w/m2s]"},
  {"SRC",     /* 198 */ "skin resevoir content [m]"},
  {"VEG",     /* 199 */ "sfc vegetation cover [%]"},
  {"VSO",     /* 200 */ "variance of subgrid scale orgography [m2]"},
  {"MX2T",    /* 201 */ "max 2m temp [K]"},
  {"MN2T",    /* 202 */ "min 2m temp [K]"},
  {"v203",    /* 203 */ "undefined"},
  {"PAW",     /* 204 */ "precip analysis weights"},
  {"RO",      /* 205 */ "runoff [m]"},
  {"v206",    /* 206 */ "undefined"},
  {"v207",    /* 207 */ "undefined"},
  {"v208",    /* 208 */ "undefined"},
  {"v209",    /* 209 */ "undefined"},
  {"v210",    /* 210 */ "undefined"},
  {"v211",    /* 211 */ "undefined"},
  {"v212",    /* 212 */ "undefined"},
  {"v213",    /* 213 */ "undefined"},
  {"v214",    /* 214 */ "undefined"},
  {"v215",    /* 215 */ "undefined"},
  {"v216",    /* 216 */ "undefined"},
  {"v217",    /* 217 */ "undefined"},
  {"v218",    /* 218 */ "undefined"},
  {"v219",    /* 219 */ "undefined"},
  {"v220",    /* 220 */ "undefined"},
  {"v221",    /* 221 */ "undefined"},
  {"v222",    /* 222 */ "undefined"},
  {"v223",    /* 223 */ "undefined"},
  {"v224",    /* 224 */ "undefined"},
  {"v225",    /* 225 */ "undefined"},
  {"v226",    /* 226 */ "undefined"},
  {"v227",    /* 227 */ "undefined"},
  {"TP",      /* 228 */ "total precip [m]"},
  {"IEWS",   /* 229 */ "instanteous sfc stress u [Nm-2]"},
  {"INSS",   /* 230 */ "instanteous sfc stress v [Nm-2]"},
  {"ISHF",   /* 231 */ "instanteous sfc sensible heat flux [Wm-2]"},
  {"IE",   /* 232 */ "instanteous sfc latent heat flux [kg/m2s]"},
  {"ASQ",    /* 233 */ "apparent sfc humidity [kg/kg]"},
  {"LSRH",  /* 234 */ "log sfc roughness"},
  {"SKT",     /* 235 */ "skin temperature [K]"},
  {"STL4",    /* 236 */ "soil temperature level 4 [K]"},
  {"SWL4",   /* 237 */ "soil wetness level 4 [m (H2O)]"},
  {"TSN",     /* 238 */ "t of snow layer [K]"},
  {"CSF",    /* 239 */ "convective snow [m]"},
  {"LSF",    /* 240 */ "large scale snow [m]"},
  {"v241",    /* 241 */ "undefined"},
  {"v242",    /* 242 */ "undefined"},
  {"FAL",   /* 243 */ "forecast albedo"},
  {"FSR",    /* 244 */ "forecast sfc roughness [m]"},
  {"FLSR",/* 245 */ "log of forecast sfc roughness"},
  {"CLWC",   /* 246 */ "Cloud liquid water content [kg/kg]"},
  {"CIWC",    /* 247 */ "Cloud ice water content [kg/kg]"},
  {"CC",    /* 248 */ "Cloud cover [0-1]"},
  {"v249",    /* 249 */ "undefined"},
  {"v250",    /* 250 */ "Ice age (0 first year, 1 multi year) [0,1]"},
  {"v251",    /* 251 */ "undefined"},
  {"v252",    /* 252 */ "undefined"},
  {"v253",    /* 253 */ "undefined"},
  {"v254",    /* 254 */ "undefined"},
  {"v255",    /* 255 */ "undefined"}
};


struct ParmTable parm_table_ecmwf_160[256] = {
  {"var0", "undefined"}, /* 0 */
  {"PRES", "Pressure [Pa]"},     /* 1 */
  {"PRMSL", "Pressure reduced to MSL [Pa]"},     /* 2 */
  {"PTEND", "Pressure tendency [Pa/s]"},     /* 3 */
  {"var4", "undefined"}, /* 4 */
  {"var5", "undefined"}, /* 5 */
  {"GP", "Geopotential [m**2/s**2]"},     /* 6 */
  {"HGT", "Geopotential height [gpm]"},     /* 7 */
  {"DIST", "Geometric height [m]"},     /* 8 */
  {"HSTDV", "Standard deviation of height [m]"},     /* 9 */
  {"HVAR", "Variance of height [m**2]"},     /* 10 */
  {"TMP", "Temperature [K]"},     /* 11 */
  {"VTMP", "Virtual temperature [K]"},     /* 12 */
  {"POT", "Potential temperature [K]"},     /* 13 */
  {"EPOT", "Pseudo-adiabatic potential temperature [K]"},     /* 14 */
  {"TMAX", "Maximum temperature [K]"},     /* 15 */
  {"TMIN", "Minimum temperature [K]"},     /* 16 */
  {"DPT", "Dew point temperature [K]"},     /* 17 */
  {"DEPR", "Dew point depression [K]"},     /* 18 */
  {"LAPR", "Lapse rate [K/m]"},     /* 19 */
  {"VISIB", "Visibility [m]"},     /* 20 */
  {"RDSP1", "Radar spectra (1) [dimensionless]"},     /* 21 */
  {"RDSP2", "Radar spectra (2) [dimensionless]"},     /* 22 */
  {"RDSP3", "Radar spectra (3) [dimensionless]"},     /* 23 */
  {"var24", "undefined"}, /* 24 */
  {"TMPA", "Temperature anomaly [K]"},     /* 25 */
  {"PRESA", "Pressure anomaly [Pa]"},     /* 26 */
  {"GPA", "Geopotential height anomaly [gpm]"},     /* 27 */
  {"WVSP1", "Wave spectra (1) [dimensionless]"},     /* 28 */
  {"WVSP2", "Wave spectra (2) [dimensionless]"},     /* 29 */
  {"WVSP3", "Wave spectra (3) [dimensionless]"},     /* 30 */
  {"WDIR", "Wind direction [degree]"},     /* 31 */
  {"WIND", "Wind speed [m/s]"},     /* 32 */
  {"UGRD", "u wind [m/s]"},     /* 33 */
  {"VGRD", "v wind [m/s]"},     /* 34 */
  {"STRM", "Stream function [m**2/s]"},     /* 35 */
  {"VPOT", "Velocity potential [m**2/s]"},     /* 36 */
  {"MNTSF", "Montgomery stream function [m**2/s**2]"},     /* 37 */
  {"SGCVV", "Sigma coord. vertical velocity [/s]"},     /* 38 */
  {"VVEL", "Pressure vertical velocity [Pa/s]"},     /* 39 */
  {"DZDT", "Geometric vertical velocity [m/s]"},     /* 40 */
  {"ABSV", "Absolute vorticity [/s]"},     /* 41 */
  {"ABSD", "Absolute divergence [/s]"},     /* 42 */
  {"RELV", "Relative vorticity [/s]"},     /* 43 */
  {"RELD", "Relative divergence [/s]"},     /* 44 */
  {"VUCSH", "Vertical u shear [/s]"},     /* 45 */
  {"VVCSH", "Vertical v shear [/s]"},     /* 46 */
  {"DIRC", "Direction of current [degree]"},     /* 47 */
  {"SPC", "Speed of current [m/s]"},     /* 48 */
  {"UOGRD", "u of current [m/s]"},     /* 49 */
  {"VOGRD", "v of current [m/s]"},     /* 50 */
  {"SPFH", "Specific humidity [kg/kg]"},     /* 51 */
  {"RH", "Relative humidity [percent]"},     /* 52 */
  {"MIXR", "Humidity mixing ratio [kg/kg]"},     /* 53 */
  {"PWAT", "Precipitable water [kg/m**2]"},     /* 54 */
  {"VAPP", "Vapor pressure [Pa]"},     /* 55 */
  {"SATD", "Saturation deficit [Pa]"},     /* 56 */
  {"EVP", "Evaporation [kg/m**2]"},     /* 57 */
  {"CICE", "Cloud Ice [kg/m**2]"},     /* 58 */
  {"PRATE", "Precipitation rate [kg/m**2/s]"},     /* 59 */
  {"TSTM", "Thunderstorm probability [percent]"},     /* 60 */
  {"APCP", "Total precipitation [kg/m**2]"},     /* 61 */
  {"NCPCP", "Large scale precipitation [kg/m**2]"},     /* 62 */
  {"ACPCP", "Convective precipitation [kg/m**2]"},     /* 63 */
  {"SRWEQ", "Snowfall rate water equivalent [kg/m**2/s]"},     /* 64 */
  {"WEASD", "Water equiv. of accum. snow depth [kg/m**2]"},     /* 65 */
  {"SNOD", "Snow depth [m]"},     /* 66 */
  {"MIXHT", "Mixed layer depth [m]"},     /* 67 */
  {"TTHDP", "Transient thermocline depth [m]"},     /* 68 */
  {"MTHD", "Main thermocline depth [m]"},     /* 69 */
  {"MTHA", "Main thermocline anomaly [m]"},     /* 70 */
  {"TCDC", "Total cloud cover [percent]"},     /* 71 */
  {"CDCON", "Convective cloud cover [percent]"},     /* 72 */
  {"LCDC", "Low level cloud cover [percent]"},     /* 73 */
  {"MCDC", "Mid level cloud cover [percent]"},     /* 74 */
  {"HCDC", "High level cloud cover [percent]"},     /* 75 */
  {"CWAT", "Cloud water [kg/m**2]"},     /* 76 */
  {"var77", "undefined"}, /* 77 */
  {"SNOC", "Convective snow [kg/m**2]"},     /* 78 */
  {"SNOL", "Large scale snow [kg/m**2]"},     /* 79 */
  {"WTMP", "Water temperature [K]"},     /* 80 */
  {"LAND", "Land-sea mask (1=land; 0=sea) [integer]"},     /* 81 */
  {"DSLM", "Deviation of sea level from mean [m]"},     /* 82 */
  {"SFCR", "Surface roughness [m]"},     /* 83 */
  {"ALBDO", "Albedo [percent]"},     /* 84 */
  {"TSOIL", "Soil temperature [K]"},     /* 85 */
  {"SOILM", "Soil moisture content [kg/m**2]"},     /* 86 */
  {"VEG", "Vegetation [percent]"},     /* 87 */
  {"SALTY", "Salinity [kg/kg]"},     /* 88 */
  {"DEN", "Density [kg/m**2]"},     /* 89 */
  {"RUNOF", "Runoff [kg/m**2]"},     /* 90 */
  {"ICEC", "Ice concentration (ice=1; no ice=0) [1/0]"},     /* 91 */
  {"ICETK", "Ice thickness [m]"},     /* 92 */
  {"DICED", "Direction of ice drift [degree]"},     /* 93 */
  {"SICED", "Speed of ice drift [m/s]"},     /* 94 */
  {"UICE", "u of ice drift [m/s]"},     /* 95 */
  {"VICE", "v of ice drift [m/s]"},     /* 96 */
  {"ICEG", "Ice growth [m]"},     /* 97 */
  {"ICED", "Ice divergence [/s]"},     /* 98 */
  {"SNOM", "Snow melt [kg/m**2]"},     /* 99 */
  {"HTSGW", "Sig height of wind waves and swell [m]"},     /* 100 */
  {"WVDIR", "Direction of wind waves [degree]"},     /* 101 */
  {"WVHGT", "Significant height of wind waves [m]"},     /* 102 */
  {"WVPER", "Mean period of wind waves [s]"},     /* 103 */
  {"SWDIR", "Direction of swell waves [degree]"},     /* 104 */
  {"SWELL", "Significant height of swell waves [m]"},     /* 105 */
  {"SWPER", "Mean period of swell waves [s]"},     /* 106 */
  {"DIRPW", "Primary wave direction [degree]"},     /* 107 */
  {"PERPW", "Primary wave mean period [s]"},     /* 108 */
  {"DIRSW", "Secondary wave direction [degree]"},     /* 109 */
  {"PERSW", "Secondary wave mean period [s]"},     /* 110 */
  {"NSWRS", "Net short wave radiation (surface) [W/m**2]"},     /* 111 */
  {"NLWRS", "Net long wave radiation (surface) [W/m**2]"},     /* 112 */
  {"NSWRT", "Net short wave radiation (top) [W/m**2]"},     /* 113 */
  {"NLWRT", "Net long wave radiation (top) [W/m**2]"},     /* 114 */
  {"LWAVR", "Long wave radiation [W/m**2]"},     /* 115 */
  {"SWAVR", "Short wave radiation [W/m**2]"},     /* 116 */
  {"GRAD", "Global radiation [W/m**2]"},     /* 117 */
  {"var118", "undefined"}, /* 118 */
  {"var119", "undefined"}, /* 119 */
  {"var120", "undefined"}, /* 120 */
  {"LHTFL", "Latent heat flux [W/m**2]"},     /* 121 */
  {"SHTFL", "Sensible heat flux [W/m**2]"},     /* 122 */
  {"BLYDP", "Boundary layer dissipation [W/m**2]"},     /* 123 */
  {"UFLX", "Zonal component of momentum flux [N/m**2]"},     /* 124 */
  {"VFLX", "Meridional component of momentum flux [N/m**2]"},     /* 125 */
  {"WMIXE", "Wind mixing energy [J]"},     /* 126 */

  /* ECMWF local table */

  {"at",      /* 127 */ "Atmospheric Tide"},
  {"bdv",     /* 128 */ "Budget Values"},
  {"zg",      /* 129 */ "Geopotential [m2 s-2]"},
  {"ta",      /* 130 */ "Temperature [K]"},
  {"ua",      /* 131 */ "U-component of Wind [ms-1]"},
  {"va",      /* 132 */ "V-component of Wind [ms-1]"},
  {"hus",     /* 133 */ "Specific Humidity [kg/kg]"},
  {"pss",     /* 134 */ "Surface Pressure [Pa]"},
  {"wa",      /* 135 */ "Vertical Velocity [Pa s-1]"},
  {"prwa",    /* 136 */ "preciptable water (vapor+drops+ice) [m]"},
  {"prw",     /* 137 */ "pecipitable water [m]"},
  {"rvort",   /* 138 */ "vorticity [s-1]"},
  {"tso1",    /* 139 */ "soil moisture level 1 [m (H20)]"},
  {"mrso1",   /* 140 */ "soil temperature level 1 [K]"},
  {"snd",     /* 141 */ "Snow Depth [m]"},
  {"prl",     /* 142 */ "Large Scale Precipitation [m]"},
  {"prc",     /* 143 */ "Convective Precipitation [m]"},
  {"prs",     /* 144 */ "Snow Fall"},
  {"bld",     /* 145 */ "Boundary Layer Dissipation [Wm-2]"},
  {"hfss",    /* 146 */ "Surface Flux of Sensible Heat [Wm-2]"},
  {"hfls",    /* 147 */ "Surface Flux of Latent Heat [Wm-2]"},
  {"v148",    /* 148 */ "undefined"},
  {"v149",    /* 149 */ "undefined"},
  {"v150",    /* 150 */ "undefined"},
  {"psl",     /* 151 */ "Mean Sea Level (MSL) Pressure [Pa]    Pa"},
  {"logpsl",  /* 152 */ "Log Surface Pressure"},
  {"v153",    /* 153 */ "undefined"},
  {"v154",    /* 154 */ "undefined"},
  {"div",     /* 155 */ "Divergence [s-1]"},
  {"zg",      /* 156 */ "Height (Geopotential) [m]"},
  {"hur",     /* 157 */ "Relative Humidity [%]"},
  {"pstn",    /* 158 */ "Tendency of Surface Pressure [Pa s-1]"},
  {"v159",    /* 159 */ "undefined"},
  {"v160",    /* 150 */ "undefined"},
  {"v161",    /* 161 */ "undefined"},
  {"v162",    /* 162 */ "undefined"},
  {"v163",    /* 163 */ "undefined"},
  {"clt",     /* 164 */ "cloud cover total [0-1]"},
  {"uas",     /* 165 */ "U-wind at 10 [ms-1]"},
  {"vas",     /* 166 */ "V-wind at 10 [ms-1]"},
  {"tas",     /* 167 */ "Temperature at 2 m [K]"},
  {"tds",     /* 168 */ "Dewpoint at 2 m [K]"},
  {"rsds",    /* 169 */ "Downward SW (sfc) [Wm-2]"},
  {"tso2",    /* 170 */ "soil temperature [K]"},
  {"mrso2",   /* 171 */ "soil wetness level 2 [m (H20)]"},
  {"lsm",     /* 172 */ "land-sea mask [(0,1)]"},
  {"sfr",     /* 173 */ "sfc roughness [m]"},
  {"albs",    /* 174 */ "Albedo [0-1]"},
  {"rlds",    /* 175 */ "Downard LW (sfc) [Wm-2]"},
  {"rss",     /* 176 */ "Net Shortwave Radiation (surface) [Wm-2]"},
  {"rls",     /* 177 */ "Net Longwave Radiation (surface) [Wm-2]"},
  {"rst",     /* 178 */ "Net Shortwave Radiation (toa) [Wm-2"},
  {"rlt",     /* 179 */ "Net Longwave Radiation (toa) [Wm-2]"},
  {"tauu",    /* 180 */ "U-component of Surface Wind Stress [Nm-2]"},
  {"tauv",    /* 181 */ "V-component of Surface Wind Stress [Nm-2]"},
  {"evs",     /* 182 */ "Evaporation [m (H2O)]"},
  {"tso3",    /* 183 */ "soil temp level 3 [m (H2O)]"},
  {"mrso3",   /* 184 */ "soil moisture level 3 [K]"},
  {"clcc",    /* 185 */ "cloud convective [0-1]"},
  {"cll",     /* 186 */ "cloud low [0-1]"},
  {"clm",     /* 187 */ "cloud mid [0-1]"},
  {"clh",     /* 188 */ "cloud high [0-1]"},
  {"v189",    /* 189 */ "undefined"},
  {"orgv",    /* 190 */ "orographic variance"},
  {"orgvew",  /* 191 */ "orographic variance e-w"},
  {"orgvns",  /* 192 */ "orographic variance n-s"},
  {"orgvnwse",/* 193 */ "orographic variance nw-se"},
  {"orgvnesw",/* 194 */ "orographic variance ne-sw"},
  {"gwsv",    /* 195 */ "gravity wave stress n-s"},
  {"gwsu",    /* 196 */ "gravity wave stress e-w"},
  {"gwd",     /* 197 */ "gravity wave diss"},
  {"src",     /* 198 */ "skin resevoir content"},
  {"sfvc",    /* 199 */ "sfc vegetation cover"},
  {"orgvsg",  /* 200 */ "orgographic variance subgrid"},
  {"tasmx",   /* 201 */ "max sfc temp"},
  {"tasmn",   /* 202 */ "min sfc temp"},
  {"v203",    /* 203 */ "undefined"},
  {"praw",    /* 204 */ "precip analysis weights"},
  {"mrro",    /* 205 */ "runoff"},
  {"cvzz",    /* 206 */ "zz variance"},
  {"cvtz",    /* 207 */ "tz covariance"},
  {"cvtt",    /* 208 */ "tt variance"},
  {"cvqz",    /* 209 */ "qz covariance"},
  {"cvqt",    /* 210 */ "qt covariance"},
  {"cvqq",    /* 211 */ "qq variance"},
  {"cvuz",    /* 212 */ "uz covariance"},
  {"cvut",    /* 213 */ "ut covariance"},
  {"cvuq",    /* 214 */ "uq covariance"},
  {"cvuu",    /* 215 */ "uu variance"},
  {"cvvz",    /* 216 */ "vz covariance"},
  {"cvvt",    /* 217 */ "vt covariance"},
  {"cvvq",    /* 218 */ "vq covariance"},
  {"cvvu",    /* 219 */ "vu covariance"},
  {"cvvv",    /* 220 */ "vv variance"},
  {"cvwz",    /* 221 */ "wz covariance"},
  {"cvwt",    /* 222 */ "wt covariance"},
  {"cvwq",    /* 223 */ "wq covariance"},
  {"cvwu",    /* 224 */ "wu covariance"},
  {"cvwv",    /* 225 */ "wv covariance"},
  {"cvww",    /* 226 */ "ww variance"},
  {"cvrr",    /* 227 */ "rh variance"},
  {"pr",      /* 228 */ "total precip"},
  {"tauui",   /* 229 */ "instanteous sfc stress u [Nm-2]"},
  {"tauvi",   /* 230 */ "instanteous sfc stress v [Nm-2]"},
  {"hfssi",   /* 231 */ "instanteous sfc sensible heat flux [Wm-2]"},
  {"hflsi",   /* 232 */ "instanteous sfc latent heat flux"},
  {"husa",    /* 233 */ "apparent sfc humidity"},
  {"logsfr",  /* 234 */ "log sfc roughness"},
  {"tgs",     /* 235 */ "skin temperature [K]"},
  {"tso4",    /* 236 */ "soil temperature level 4 [K]"},
  {"mrso4",   /* 237 */ "soil wetness level 4 [m (H2O)]"},
  {"tgs",     /* 238 */ "t of snow layer [K]"},
  {"prsc",    /* 239 */ "convective snow [m]"},
  {"prsl",    /* 240 */ "large scale snow [m]"},
  {"cllw",    /* 241 */ "cloud liquid water"},
  {"clct",    /* 242 */ "total cloud cover"},
  {"albsf",   /* 243 */ "forecast albedo"},
  {"sfrf",    /* 244 */ "forecast sfc roughness"},
  {"logsfcrf",/* 245 */ "log offorecast sfc roughness"},
  {"wspds",   /* 246 */ "10 m wind speed"},
  {"taum",    /* 247 */ "magnitude of momentum flux"},
  {"v248",    /* 248 */ "undefined"},
  {"gwmf",    /* 249 */ "gravity wave drag momentum flux [Nm-2]"},
  {"v250",    /* 250 */ "undefined"},
  {"v251",    /* 251 */ "undefined"},
  {"v252",    /* 252 */ "undefined"},
  {"v253",    /* 253 */ "undefined"},
  {"v254",    /* 254 */ "undefined"},
  {"v255",    /* 255 */ "undefined"}
};

/* parameter table for ocean modeling branch (OMB) of NCEP */
/* center = 7, subcenter = EMC, parameter table = 128 */

struct ParmTable parm_table_omb[256] = {
 {"var0", "Reserved"},
 {"var1", "Reserved"},
 {"GHz6", "6.6 GHz - K"},
 {"GHz10", "10.7 GHz - K"},
 {"GHz18", "18.0 GHz - K"},
 {"GHz19V", "SSMI 19 GHz, Vertical Polarization - K"},
 {"GHz19H", "SSMI 19 GHz, Horizontal Polarization - K"},
 {"GHz21", "21.0 GHz - K"},
 {"GHz22V", "SSMI 22 GHz, Vertical Polarization - K"},
 {"GHz37V", "SSMI 37 GHz, Vertical Polarization - K"},
 {"GHz37H", "SSMI 37 GHz, Horizontal Polarization - K"},
 {"MSU1", "MSU Ch 1 - 50.30 GHz - K"},
 {"MSU2", "MSU Ch 2 - 53.74 GHz - K"},
 {"MSU3", "MSU Ch 3 - 54.96 GHz - K"},
 {"MSU4", "MSU Ch 4 - 57.95 GHz - K"},
 {"GHz85V", "SSMI 85 GHz, Vertical Polarization - K"},
 {"GHz85H", "SSMI 85 GHz, Horizontal Polarization - K"},
 {"GHz91", "91.65 GHz - K"},
 {"GHz150", "150 GHz - K"},
 {"GHz183pm7", "183 +- 7 GHz - K"},
 {"GHz183pm3", "183 +- 3 GHz - K"},
 {"GHz183pm1", "183 +- 1 GHz - K"},
 {"SSMT1C1", "SSM/T1 - ch 1 - K"},
 {"SSMT1C2", "SSM/T1 - ch 2 - K"},
 {"SSMT1C3", "SSM/T1 - ch 3 - K"},
 {"SSMT1C4", "SSM/T1 - ch 4 - K"},
 {"SSMT1C5", "SSM/T1 - ch 5 - K"},
 {"SSMT1C6", "SSM/T1 - ch 6 - K"},
 {"SSMT1C7", "SSM/T1 - ch 7 - K"},
 {"var29", "Reserved"},
 {"var30", "Reserved"},
 {"var31", "Reserved"},
 {"var32", "Reserved"},
 {"var33", "Reserved"},
 {"var34", "Reserved"},
 {"var35", "Reserved"},
 {"var36", "Reserved"},
 {"var37", "Reserved"},
 {"var38", "Reserved"},
 {"var39", "Reserved"},
 {"var40", "Reserved"},
 {"var41", "Reserved"},
 {"var42", "Reserved"},
 {"var43", "Reserved"},
 {"var44", "Reserved"},
 {"var45", "Reserved"},
 {"var46", "Reserved"},
 {"var47", "Reserved"},
 {"var48", "Reserved"},
 {"var49", "Reserved"},
 {"var50", "Reserved"},
 {"var51", "Reserved"},
 {"var52", "Reserved"},
 {"var53", "Reserved"},
 {"var54", "Reserved"},
 {"var55", "Reserved"},
 {"var56", "Reserved"},
 {"var57", "Reserved"},
 {"var58", "Reserved"},
 {"var59", "Reserved"},
 {"MI14.95", "HIRS/2 ch 1 - 14.95 micron - K"},
 {"MI14.71", "HIRS/2, GOES 14.71 micron - K"},
 {"MI14.49", "HIRS/2 ch 3 - 14.49 micron - K"},
 {"MI14.37", "GOES I-M - 14.37 micron - K"},
 {"MI14.22", "HIRS/2 ch 4 - 14.22 micron - K"},
 {"MI14.06", "GOES I-M - 14.06 micron - K"},
 {"MI13.97", "HIRS/2 ch 5 - 13.97 micron - K"},
 {"MI13.64", "HIRS/2, GOES 13.64 micron - K"},
 {"MI13.37", "GOES I-M - 13.37 micron - K"},
 {"MI13.35", "HIRS/2 ch 7 - 13.35 micron - K"},
 {"MI12.66", "GOES I-M - 12.66 micron - K"},
 {"MI12.02", "GOES I-M - 12.02 micron - K"},
 {"MI12.00", "AVHRR ch 5 - 12.0 micron - K"},
 {"MI11.11", "HIRS/2 ch 8 - 11.11 micron - K"},
 {"MI11.03", "GOES I-M - 11.03 micron - K"},
 {"MI10.80", "AVHRR ch 4 - 10.8 micron - K"},
 {"MI9.71", "HIRS/2, GOES - 9.71 micron - K"},
 {"var77", "Reserved"},
 {"var78", "Reserved"},
 {"var79", "Reserved"},
 {"MI8.16", "HIRS/2 ch 10 - 8.16 micron - K"},
 {"MI7.43", "GOES I-M - 7.43 micron - K"},
 {"MI7.33", "HIRS/2 ch 11 - 7.33 micron - K"},
 {"MI7.02", "GOES I-M - 7.02 micron - K"},
 {"MI6.72", "HIRS/2 ch 12 - 6.72 micron - K"},
 {"MI6.51", "GOES I-M - 6.51 micron - K"},
 {"MI4.57", "HIRS/2, GOES - 4.57 micron - K"},
 {"MI4.52", "HIRS/2, GOES - 4.52 micron - K"},
 {"MI4.46", "HIRS/2 ch 15 - 4.46 micron - K"},
 {"MI4.45", "GOES I-M - 4.45 micron - K"},
 {"MI4.40", "HIRS/2 ch 16 - 4.40 micron - K"},
 {"MI4.24", "HIRS/2 ch 17 - 4.24 micron - K"},
 {"MI4.13", "GOES I-M - 4.13 micron - K"},
 {"MI4.00", "HIRS/2 ch 18 - 4.00 micron - K"},
 {"MI8.16", "GOES I-M - 3.98 micron - K"},
 {"MI8.16", "HIRS/2 Window - 3.76 micron - K"},
 {"MI8.16", "AVHRR, GOES - 3.74 micron - K"},
 {"var97", "Reserved"},
 {"var98", "Reserved"},
 {"var99", "Reserved"},
 {"MI0.91", "AVHRR ch 2 - 0.91 micron - K"},
 {"MI0.696", "GOES I-M - 0.696 micron - K"},
 {"MI0.69", "HIRS/2 Vis - 0.69 micron - K"},
 {"MI0.63", "AVHRR ch 1 - 0.63 micron - K"},
 {"var104", "Reserved"},
 {"var105", "Reserved"},
 {"var106", "Reserved"},
 {"var107", "Reserved"},
 {"var108", "Reserved"},
 {"var109", "Reserved"},
 {"var110", "Reserved"},
 {"var111", "Reserved"},
 {"var112", "Reserved"},
 {"var113", "Reserved"},
 {"var114", "Reserved"},
 {"var115", "Reserved"},
 {"var116", "Reserved"},
 {"var117", "Reserved"},
 {"var118", "Reserved"},
 {"var119", "Reserved"},
 {"var120", "Reserved"},
 {"var121", "Reserved"},
 {"var122", "Reserved"},
 {"var123", "Reserved"},
 {"var124", "Reserved"},
 {"var125", "Reserved"},
 {"var126", "Reserved"},
 {"var127", "Reserved"},
 {"AVDEPTH", "Ocean depth - mean - m"},
 {"DEPTH", "Ocean depth - instantaneous - m"},
 {"ELEV", "Ocean surface elevation relative to geoid - m"},
 {"MXEL24", "Max ocean surface elevation in last 24 hours - m"},
 {"MNEL24", "Min ocean surface elevation in last 24 hours - m"},
 {"var133", "Reserved"},
 {"var134", "Reserved"},
 {"O2", "Oxygen -Mol/kg"},
 {"PO4", "PO4 - Mol/kg"},
 {"NO3", "NO3 - Mol/kg"},
 {"SiO4", "SiO4 - Mol/kg"},
 {"CO2aq", "CO2 (aq) - Mol/kg"},
 {"HCO3", "HCO3 - - Mol/kg"},
 {"CO3", "CO3 -- - Mol/kg"},
 {"TCO2", "TCO2 - Mol/kg"},
 {"TALK", "TALK - Mol/kg"},
 {"var144", "Reserved"},
 {"var145", "Reserved"},
 {"S11", "S11 - 1,1 component of ice stress tensor"},
 {"S12", "S12 - 1,2 component of ice stress tensor"},
 {"S22", "S22 - 2,2 component of ice stress tensor"},
 {"INV1", "T1 - First invariant of stress tensor"},
 {"INV2", "T2 - Second invariant of stress tensor"},
 {"var151", "Reserved"},
 {"var152", "Reserved"},
 {"var153", "Reserved"},
 {"var154", "Reserved"},
 {"WVRGH", "Wave Roughness"},
 {"WVSTRS", "Wave Stresses"},
 {"WHITE", "Whitecap coverage"},
 {"SWDIRWID", "Swell direction width"},
 {"SWFREWID", "Swell frequency width"},
 {"WVAGE", "Wave age"},
 {"PWVAGE", "Physical Wave age"},
 {"var162", "Reserved"},
 {"var163", "Reserved"},
 {"var164", "Reserved"},
 {"LTURB", "Master length scale (turbulence) - m"},
 {"var166", "Reserved"},
 {"var167", "Reserved"},
 {"var168", "Reserved"},
 {"var169", "Reserved"},
 {"AIHFLX", "Net Air-Ice heat flux - W/m^2"},
 {"AOHFLX", "Net Air-Ocean heat flux - W/m^2"},
 {"IOHFLX", "Net Ice-Ocean heat flux - W/m^2"},
 {"IOSFLX", "Net Ice-Ocean salt flux - kg/s"},
 {"var174", "Reserved"},
 {"OMLT", "Ocean Mixed Layer Temperature - K"},
 {"OMLS", "Ocean Mixed Layer Salinity - kg/kg"},
 {"var177", "Reserved"},
 {"var178", "Reserved"},
 {"var179", "Reserved"},
 {"var180", "Reserved"},
 {"var181", "Reserved"},
 {"var182", "Reserved"},
 {"var183", "Reserved"},
 {"var184", "Reserved"},
 {"var185", "Reserved"},
 {"var186", "Reserved"},
 {"var187", "Reserved"},
 {"var188", "Reserved"},
 {"var189", "Reserved"},
 {"var190", "Reserved"},
 {"var191", "Reserved"},
 {"var192", "Reserved"},
 {"var193", "Reserved"},
 {"var194", "Reserved"},
 {"var195", "Reserved"},
 {"var196", "Reserved"},
 {"var197", "Reserved"},
 {"var198", "Reserved"},
 {"var199", "Reserved"},
 {"var200", "Reserved"},
 {"var201", "Reserved"},
 {"var202", "Reserved"},
 {"var203", "Reserved"},
 {"var204", "Reserved"},
 {"var205", "Reserved"},
 {"var206", "Reserved"},
 {"var207", "Reserved"},
 {"var208", "Reserved"},
 {"var209", "Reserved"},
 {"var210", "Reserved"},
 {"var211", "Reserved"},
 {"var212", "Reserved"},
 {"var213", "Reserved"},
 {"var214", "Reserved"},
 {"var215", "Reserved"},
 {"var216", "Reserved"},
 {"var217", "Reserved"},
 {"var218", "Reserved"},
 {"var219", "Reserved"},
 {"var220", "Reserved"},
 {"var221", "Reserved"},
 {"var222", "Reserved"},
 {"var223", "Reserved"},
 {"var224", "Reserved"},
 {"var225", "Reserved"},
 {"var226", "Reserved"},
 {"var227", "Reserved"},
 {"var228", "Reserved"},
 {"var229", "Reserved"},
 {"var230", "Reserved"},
 {"var231", "Reserved"},
 {"var232", "Reserved"},
 {"var233", "Reserved"},
 {"var234", "Reserved"},
 {"var235", "Reserved"},
 {"var236", "Reserved"},
 {"var237", "Reserved"},
 {"var238", "Reserved"},
 {"var239", "Reserved"},
 {"var240", "Reserved"},
 {"var241", "Reserved"},
 {"var242", "Reserved"},
 {"var243", "Reserved"},
 {"var244", "Reserved"},
 {"var245", "Reserved"},
 {"var246", "Reserved"},
 {"var247", "Reserved"},
 {"var248", "Reserved"},
 {"var249", "Reserved"},
 {"var250", "Reserved"},
 {"var251", "Reserved"},
 {"var252", "Reserved"},
 {"var253", "Reserved"},
 {"var254", "Reserved"},
 {"var255", "Reserved"}
};


/*
 * EC_ext	v1.0 wesley ebisuzaki
 *
 * prints something readable from the EC stream parameter
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


/*
 * get grid size from GDS
 *
 */

int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny) {

    int i, ix, iy;
    long int isum;

    *nx = ix = GDS_LatLon_nx(gds);
    *ny = iy = GDS_LatLon_ny(gds);
    *nxny = ix * iy;

    if (GDS_LatLon(gds) | GDS_Mercator(gds) | GDS_Gnomonic(gds)) return 0;

    if (GDS_Gaussian(gds)) {
	if (ix == 65535) {
	    *nx = -1;
	    /* reduced grid */
	    isum = 0;
	    for (i = 0; i < iy / 2; i++) {
		isum += gds[32+i*2]*256 + gds[33+i*2];
	    }
	    isum += isum;
	    if (iy % 2) {
		i = iy / 2 + 1;
		isum += gds[32+i*2]*256 + gds[33+i*2];
	    }
	    *nxny = isum;
	}
	return 0;
    }
    return 0;
}

#define NCOL 15
void GDS_prt_thin_lon(unsigned char *gds) {
    int iy, i, col;

    iy = GDS_LatLon_ny(gds);
    iy = (iy + 1) / 2;

    for (col = i = 0; i < iy; i++) {
	if (col == 0) printf("   ");
	printf("%5d", (gds[32+i*2] << 8) + gds[33+i*2]);
	col++;
	if (col == NCOL) {
	    col = 0;
	    printf("\n");
	}
    }
    if (col != 0) printf("\n");
}

#define START -1

static int user_center = 0, user_subcenter = 0, user_ptable = 0;
static enum {filled, not_found, not_checked, no_file, init} status = init;

struct ParmTable parm_table_user[256];

/*
 * sets up user parameter table
 */

int setup_user_table(int center, int subcenter, int ptable) {

    int i, j, c0, c1, c2;
    FILE *input;
    char *filename, line[300];

    if (status == init) {
	for (i = 0; i < 256; i++) {
	    parm_table_user[i].name = parm_table_user[i].comment = NULL;
	}
	status = not_checked;
    }

    if (status == no_file) return 0;

    if ((user_center == -1 || center == user_center) &&
	    (user_subcenter == -1 || subcenter == user_subcenter) &&
	    (user_ptable == -1 || ptable == user_ptable)) {

	if (status == filled) return 1;
	if (status == not_found) return 0;
    }

    /* open gribtab file */

    filename = getenv("GRIBTAB");
    if (filename == NULL) filename = getenv("gribtab");
    if (filename == NULL) filename = "gribtab";

    if ((input = fopen(filename,"r")) == NULL) {
        status = no_file;
        return 0;
    }

    user_center = center;
    user_subcenter = subcenter;
    user_ptable = ptable;

    /* scan for center & subcenter and ptable */
    for (;;) {
        if (fgets(line, 299, input) == NULL) {
	    status = not_found;
            return 0;
        }
	if (atoi(line) != START) continue;
	i = sscanf(line,"%d:%d:%d:%d", &j, &center, &subcenter, &ptable);
        if (i != 4) {
	    fprintf(stderr,"illegal gribcap center/subcenter/ptable line: %s\n", line);
            continue;
        }
	if ((center == -1 || center == user_center) &&
	    (subcenter == -1 || subcenter == user_subcenter) &&
	    (ptable == -1 || ptable == user_ptable)) break;
    }

    user_center = center;
    user_subcenter = subcenter;
    user_ptable = ptable;

    /* free any used memory */
    if (parm_table_user[i].name != NULL) {
        for (i = 0; i < 256; i++) {
	    free(parm_table_user[i].name);
            free(parm_table_user[i].comment);
        }
    }

    /* read definitions */

    for (;;) {
        if (fgets(line, 299, input) == NULL) break;
	if ((i = atoi(line)) == START) break;
	line[299] = 0;

	/* find the colons and end-of-line */
	for (c0 = 0; line[c0] != ':' && line[c0] != 0; c0++) ;
        if (line[c0] == 0) {
	    fprintf(stderr,"illegal gribcap line:%s\n", line);
	    continue;
	}
	for (c1 = c0 + 1; line[c1] != ':' && line[c1] != 0; c1++) ;
	c2 = strlen(line);
        if (line[c2-1] == '\n') line[--c2] = '\0';
        if (c2 <= c1) {
	    fprintf(stderr,"illegal gribcap line:%s\n", line);
	    continue;
	}
	line[c0] = 0;
	line[c1] = 0;

	parm_table_user[i].name = (char *) malloc(c1 - c0);
	parm_table_user[i].comment = (char *) malloc(c2 - c1);
	strcpy(parm_table_user[i].name, line+c0+1);
	strcpy(parm_table_user[i].comment, line+c1+1);
    }

    /* now to fill in undefined blanks */
    for (i = 0; i < 255; i++) {
	if (parm_table_user[i].name == NULL) {
	    parm_table_user[i].name = (char *) malloc(7);
	    sprintf(parm_table_user[i].name, "var%d", i);
	    parm_table_user[i].comment = (char *) malloc(strlen("undefined")+1);
	    strcpy(parm_table_user[i].comment, "undefined");
        }
    }
    status = filled;
    return 1;
}
