#ifndef DRIVER_GAGMAP
#define WHERE extern
#else
#define WHERE
#endif

WHERE FILE *gfile;
WHERE FILE *mfile;

/* Following structure holds all the unpacked header info from
   a grib record. */

struct grhdr {
  int vers;
  int len;
  int pdslen,gdslen,bmslen,bdslen;
  int id;
  int gdsflg,bmsflg;
  int parm;
  int ltyp;
  int level;
  int l1,l2;
  struct dt dtim;
  struct dt btim;
  int ftu,p1,p2,tri;
  int fcstu,fcstt;
  int cent;
  float dsf;
  int gtyp,gicnt,gjcnt,gsf1,gsf2,gsf3;
  int bnumr;
  int bpos;
  int iflg;
  float bsf;
  float ref;
  int bnum;
  int dpos;
};


/* ---------------- global variables ------------------- */

WHERE int fpos;           /* File pointer into GRIB file */
WHERE int verb;           /* Verbose option */
WHERE int quiet;          /* quite option */
WHERE int ver;            /* version */
WHERE int diag;           /* Verbose option */
WHERE int irec;
WHERE int scanflg;        /* general scan between GRIB records ASSUMED */
WHERE int scaneof;        /* option to ignore failure to find data at end of file */
WHERE int scanEOF;        /* option to ignore failure to find data at end of file */
WHERE int scanlim;        /* the default # of max bytes between records */
WHERE int notau;          /* force time to be base time */
WHERE int tauflg;         /* search for a fixed tau in filling the 4-D volume */
WHERE int tauoff;         /* the fixed tau in h */
WHERE int tau0;           /* set the base dtg for tau search */
WHERE int forceok;        /* set the base dtg for tau search */
WHERE int mpiflg;         /* Artificial initial date/time same as tau0!!*/
WHERE int write_map;      /* write out the map (testing only) */
WHERE int update;         /* update mode for templated files for NCEP CPC */
WHERE struct dt btimdd;   /* initial base time from dd file */

WHERE int nrec;           /* Number of records per grid */
WHERE int gtype[16];      /* Grid types for this grid set */
WHERE struct gafile *pfi;

WHERE struct gaindx *pindx;
WHERE struct dt dtim, dtimi;
WHERE int cnt,rc,i,flg,iarg,tcur,told;
WHERE char cmd[256];
WHERE char rec[512], *ch, *ifile;

WHERE int len, skip;
WHERE struct grhdr ghdr;


/* ---------------- prototypes ------------------- */

extern int gribhdr(struct grhdr *);
extern int gribrec(struct grhdr *, struct gafile *, struct gaindx *);
extern void gribfill (int, int, int, struct grhdr *, struct gaindx *);
extern void gribpr (struct grhdr *);
extern int gribmap (void) ;

