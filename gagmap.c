#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

#include "grads.h"
#include "gagmap.h"

/*mf 961205 --- expose Mike Fiorino's global struct to these routines for 365 day calandars mf*/
extern struct gamfcmn mfcmn;
/*mf 961205 --- expose Mike Fiorino's global struct to these routines for 365 day calandars mf*/

/*  Routine to scan a grib file and output a grib index (map) file. */

/*---- 960829 Modification to run outside of gribmap for LATS

   Mike Fiorino to run outside gribmap 

----*/

/*  Values output into the grib map file:

    Header:

    hipnt info:  0 - version number (1)
                 1 - number of times in file
                 2 - number of records per time
                 3 - Grid type
                   255 - user defined grid.  descriptor
                         describes grid exactly; one record
                         per grid.
                    29 - Predefined grid set 29 and 30.
                         Two records per grid.

    hfpnt info:  None

    Info:

    intpnt info (for each mapped grib record) :
                 0 - position of start of data in file
                 1 - position of start of bit map in file
                 2 - number of bits per data element

    fltpnt info :
                 0 - decimal scale factor for this record
                 1 - binary scale factor
                 2 - reference value

*/


struct dt rtime;   /* time based on input map */
struct dt ftime;   /* time based on dd file */

struct gaindx *rindx;


/* internal prototypes */

int wtgmap(void) ;
void putint(int, unsigned char *,int *) ;

/* global to this file only */

int nelements=3;
char mrec[512];
static int flen;
/*--------------------------------------- 

  main routine

---------------------------------------*/

int gribmap (void) {


/*---
  declare internal vars
---*/

long filesize;
unsigned char vermap;
int idum;
float fdum;
unsigned char urec[4];
int nt,i,j;
int update_map=0;

float map_grvals[8];
int rt,rtoff,toff,nxytime;
int didmatch=0,rcgr;

/*---
  initialize this global variable here... 
---*/

mfile=NULL;

/*--- 
  open the dd file
---*/

  if (ifile==NULL) {
    printf ("\n");
    cnt = nxtcmd (cmd,"Enter name of Data Descriptor file: ");
    if (cnt==0) return(1);
    getwrd(rec,cmd,250);
    ifile = rec;
  }

  pfi = getpfi();
  if (pfi==NULL) {
   printf ("Memory allocation error \n");
   return(1);
  }

  pindx = (struct gaindx *)malloc(sizeof(struct gaindx));
  if (pindx==NULL) {
    printf ("Memory allocation error for pindx\n");
    return(1);
  }

  rindx = (struct gaindx *)malloc(sizeof(struct gaindx));
  if (rindx==NULL) {
    printf ("Memory allocation error for rindx\n");
    return(1);
  }

/* ---------------- open the dd file  ------------------- */


  rc = gaddes (ifile, pfi, 0);
  if (rc) {
    if(! quiet) printf ("File name is:  %s\n",ifile);
    return(1);
  }
  if (pfi->idxflg!=1) {
    printf ("Data Descriptor file is not for GRIB data\n");
    return(1);
  }

  /* save the initial time from the data descriptor file for the tau0 option and the map file*/

  btimdd.yr = *(pfi->abvals[3]);
  btimdd.mo = *(pfi->abvals[3]+1);
  btimdd.dy = *(pfi->abvals[3]+2);
  btimdd.hr = *(pfi->abvals[3]+3);
  btimdd.mn = *(pfi->abvals[3]+4);

  if(diag) printf(" btim %i %i %i %i %i\n",btimdd.yr,btimdd.mo,btimdd.dy,btimdd.hr,btimdd.mn);


  /* Set up for this grid type */

  if (pfi->grbgrd<-900 || pfi->grbgrd==255) {
    nrec = 1;
    gtype[0] = 255;
  } else if (pfi->grbgrd>-1 && pfi->ppflag) {
    gtype[0] = pfi->grbgrd;
    nrec=1;
  } else if (pfi->grbgrd==29) {
    nrec = 2;
    gtype[0] = 29;
    gtype[1] = 30;
    if (pfi->dnum[0]!=144 || pfi->dnum[1]!=73 ||
        pfi->linear[0]!=1 || pfi->linear[1]!=1 ||
        *(pfi->grvals[0])!= 2.5 || *(pfi->grvals[0]+1) != -2.5 ||
        *(pfi->grvals[1])!= 2.5 || *(pfi->grvals[1]+1) != -92.5 ) {
      printf ("Error in grid specification for GRIB grid type 29/30.\n");
      printf ("  Grid scaling must indicate a 2.5 x 2.5 grid\n");
      printf ("  Grid size must be 144 x 73\n");
      printf ("  Grid must go from 0 to 357.5 and -90 to 90\n");
      return(1);
    }
  } else {
    nrec = 1;
    gtype[0] = pfi->grbgrd;
  }

  /* Set index stuff */

  pindx->type = ver;  /* gribmap version */
  pindx->hinum = 4;
  pindx->hfnum = 0;
  pindx->intnum = nrec * nelements * pfi->trecs * pfi->dnum[3];
  pindx->fltnum = nrec * nelements * pfi->trecs * pfi->dnum[3];
  pindx->hipnt = (int *)malloc(sizeof(int)*pindx->hinum);
  pindx->intpnt = (int *)malloc(sizeof(int)*pindx->intnum);
  pindx->fltpnt = (float *)malloc(sizeof(float)*pindx->fltnum);

  if (pindx->hipnt==NULL || pindx->intpnt==NULL || pindx->fltpnt==NULL) {
    printf ("memory allocation error\n");
    return(1);
  }
for (i=0; i<pindx->intnum; i++) *(pindx->intpnt+i) = -999;
for (i=0; i<pindx->fltnum; i++) *(pindx->fltpnt+i) = -999; /* try initializing */
*(pindx->hipnt) = ver;
*(pindx->hipnt+1) = pfi->dnum[3];
*(pindx->hipnt+2) = pfi->trecs;
*(pindx->hipnt+3) = pfi->grbgrd;
if (pfi->grbgrd<-900) *(pindx->hipnt+3) = 255;

if(! quiet) printf ("\ngribmap:  Scanning binary GRIB file(s):\n");


/* If there are multiple files in a time series, we need to
   loop through them based on the info in pfi */

tcur = 0;
gfile = NULL;

if(pfi->tmplat && update ) {

  if(!quiet) {
    printf("UUUUU: Updating the gribmap for the TEMPLATED GRIB data described by the:\n");
    printf("UUUUU:\n");
    printf("UUUUU: %s\n",pfi->dnam);
    printf("UUUUU:\n");
    printf("UUUUU: GrADS data descriptor (.ctl) file....\n\n");
  }

  if (nrec != 1) {
    printf ("Attempting to update a special case GRIB file nrec = %d\n",nrec);
    printf ("No can do, bye \n");
    return(42);
  }

  mfile = fopen(pfi->mnam,"rb+");

  if (mfile==NULL) {

    printf ("Could not open index file: %s for updating\n",pfi->mnam);
    printf ("Create the map........\n");
    update=0;

  } else {

    fseek(mfile,0L,SEEK_END);
    filesize=ftell(mfile);
  
    if(filesize>0) {

      fseek(mfile,1,0);
      fread(&vermap,sizeof(unsigned char),1,mfile);
      if(diag) printf(" vermap =%d\n",vermap);

      if(vermap == 2) {

	fseek(mfile,2,0);
	fread(mrec,sizeof(unsigned char),4,mfile);
	rindx->hinum=gagby(mrec,0,4);
	if(diag) printf(" hinum = %d\n",rindx->hinum);
	
	fread(mrec,sizeof(unsigned char),4,mfile);
	rindx->hfnum=gagby(mrec,0,4);
	if(diag) printf(" hfnum = %d\n",rindx->hfnum);

	fread(mrec,sizeof(unsigned char),4,mfile);
	rindx->intnum=gagby(mrec,0,4);
	if(diag) printf(" intnum = %d\n",rindx->intnum);

	fread(mrec,sizeof(unsigned char),4,mfile);
	rindx->fltnum=gagby(mrec,0,4);
	if(diag) printf(" fltnum = %d\n",rindx->fltnum);

	rindx->hipnt = NULL;
	rindx->hfpnt = NULL;
	rindx->intpnt = NULL;
	rindx->fltpnt = NULL;

/* skip the begining time struct info */

	fread(mrec,sizeof(unsigned char),7,mfile);

	if (rindx->hinum>0) {
	  rindx->hipnt = (int *)malloc(sizeof(int)*rindx->hinum);
	  if (rindx->hipnt==NULL) goto err8;

	  for(i=0;i<rindx->hinum;i++) { 
	    fread(mrec,sizeof(unsigned char),4,mfile);
	    idum=gagby(mrec,0,4);
	    if(gagbb(mrec,0,1)) idum=-idum;
	    *(rindx->hipnt+i)=idum;
	    if(diag) printf(" 1 i = %d rindx->hipnt = %d\n",i,idum);
	  }

	}

	if (rindx->hfnum>0) {
	  rindx->hfpnt = (float *)malloc(sizeof(float)*rindx->hfnum);
	  if (rindx->hfpnt==NULL) goto err8;
	  fread (rindx->hfpnt,sizeof(float),rindx->hfnum,mfile);
	}

	if (rindx->intnum>0) {
	  rindx->intpnt = (int *)malloc(sizeof(int)*rindx->intnum);
	  if (rindx->intpnt==NULL) goto err8;
	  for(i=0;i<rindx->intnum;i++) { 
	    fread(mrec,sizeof(unsigned char),4,mfile);
	    idum=gagby(mrec,0,4);
	    if(gagbb(mrec,0,1)) idum=-gagbb(mrec,1,31);
	    *(rindx->intpnt+i)=idum;
	    if(diag) printf(" 2 i = %d rindx->intpnt = %d\n",i,idum);
	  }
	  
	}
	
	if (rindx->fltnum>0) {
	  rindx->fltpnt = (float *)malloc(sizeof(float)*rindx->fltnum);
	  if (rindx->fltpnt==NULL) goto err8;
	  
	  for(i=0;i<rindx->fltnum;i++) { 
	    fread(urec,sizeof(unsigned char),4,mfile);
	    fdum=ibm2flt(urec);
	    *(rindx->fltpnt+i)=fdum;
	    if(diag) printf(" 3 i = %d rindx->fltpnt = %g\n",i,fdum);
	  }

	}

/* read in the grid to absolute time factors at the end of the regular map info */

	for(i=0;i<8;i++) { 
	  fread(urec,sizeof(unsigned char),4,mfile);
	  fdum=ibm2flt(urec);
	  *(map_grvals+i)=fdum;
	  if(diag) printf(" grvals i = %d map_grvals = %g\n",i,fdum);
	}

      }
  

      nt=*(rindx->hipnt+1);
      nxytime=pfi->trecs*nrec*nelements;

      for(i=1;i<=nt;i++) {
	gr2t(map_grvals,(float)i,&rtime);
	rt=(int)(t2gr(map_grvals,&rtime)+0.01);
	rtoff=(rt-1)*nxytime;
	if(diag) printf("rrrr %d %d %d %d %d rt = %d rtoff = %d\n",
			rtime.yr,rtime.mo,rtime.dy,rtime.hr,rtime.mn,rt,rtoff);
      }
    

    } else {
      
      printf("WARNING:  Map file for updating has 0 size\n");
      printf("          will go ahead and do the map...\n\n");
      update=0;
      
    }
 
  }

/*--- 
  close the map file now that we've gotten all the information we need; including the map
---*/

  fclose(mfile);
  mfile=NULL;

}


while (1) {

  if (pfi->tmplat) {

    if (gfile) fclose(gfile);

    told = tcur;

    if (tcur==0) {

      tcur=1;

      if(update) {

	tcur=1;
	told=tcur;

	while (told==*(pfi->fnums+tcur-1) && tcur<=pfi->dnum[3]) {

	  gr2t(pfi->grvals[3],(float)tcur,&dtim);
	  rt=t2gr(map_grvals,&dtim);
	  rt=(int)(t2gr(map_grvals,&dtim)+0.01);
	  rtoff=(rt-1)*nxytime;
	  toff=(tcur-1)*nxytime;

/*---
  check if a match for first time in the file 
---*/

	  update_map=0;
	  if( rt > 0 && rt <= nt ) update_map=1;

	  if(diag) {
	    printf("\nuuu0 told= %d tcur= %d dtim:  %d %d %d %d %d\n",
		   told,tcur,dtim.yr,dtim.mo,dtim.dy,dtim.hr,dtim.mn);
	    printf("uuu0   rt = %d rtoff = %d toff = %d \n",rt,rtoff,toff);
	    printf("uuu0 told= %d tcur= %d update_map = %d\n",told,tcur,update_map);
	  }

/*---------------------
  check if update not available in MIDDLE of a time chunk (told -> tcur)
------------------*/

	  if(update_map==0 && (tcur != told) ) {
	    if(diag) printf("\nmmmmmmm 000  check told = %d tcur=%d\n\n",told,tcur);
	    tcur=told;
	    update=0;
	    break;
	  }

	  if(update_map) {
	    for(j=0;j<nxytime;j++) {
	      *(pindx->intpnt+toff+j)=*(rindx->intpnt+rtoff+j);
	      *(pindx->fltpnt+toff+j)=*(rindx->fltpnt+rtoff+j);
	    }
	  }

	  tcur++;

	}

	tcur=1;
      }

    } else {

      while (told==*(pfi->fnums+tcur-1) && tcur<=pfi->dnum[3]) {

	if(update) {

	  gr2t(pfi->grvals[3],(float)tcur,&dtim);
	  rt=(int)(t2gr(map_grvals,&dtim)+0.01);
	  rtoff=(rt-1)*nxytime;
	  toff=(tcur-1)*nxytime;

	  update_map=0;
	  if( rt > 0 && rt <= nt ) update_map=1;

	  if(diag) {
	    printf("\nuuu1 told= %d tcur= %d dtim:  %d %d %d %d %d\n",told,tcur,dtim.yr,dtim.mo,dtim.dy,dtim.hr,dtim.mn);
	    printf("uuu1   rt = %d  rtoff = %d toff = %d\n",rt,rtoff,toff);
	    printf("uuu1 update_map = %d\n",update_map);
	  }
	  
/*---------------------
  check if update not available in MIDDLE of a time chunk (told -> tcur)
------------------*/

	  if(update_map==0 && (tcur != told) ) {
	    if(diag) printf("\nmmmmmmm check told = %d tcur=%d\n\n",told,tcur);
	    tcur=told;
	    update=0;
	    break;
	  }

	  if(update_map) {
	    for(j=0;j<nxytime;j++) {
	      *(pindx->intpnt+toff+j)=*(rindx->intpnt+rtoff+j);
	      *(pindx->fltpnt+toff+j)=*(rindx->fltpnt+rtoff+j);
	    }
	  }

	}

	tcur++;

      }

/*---
  tcur has been updated; check if we are at the end of the file...
---*/

      if(update) {
	gr2t(pfi->grvals[3],(float)tcur,&dtim);
	rt=(int)(t2gr(map_grvals,&dtim)+0.01);
	update_map=0;
	if( rt > 0 && rt <= nt ) update_map=1;
      }

    }

    if (tcur>pfi->dnum[3]) break;
    if (tcur != *(pfi->fnums+tcur-1)) {
      printf ("Logic error in file template structure\n");
      return(1);
    }
    
    gr2t(pfi->grvals[3], (float)tcur, &dtim);
    gr2t(pfi->grvals[3], 1.0, &dtimi);

    if(diag) printf("11111 Before open test:  tcur= %d update= %d update_map=%d\n",
		    tcur,update,update_map);

    if(update_map==0) {
      ch = gafndt(pfi->name, &dtim, &dtimi, pfi->abvals[3]);
      if (ch==NULL) {
	printf ("Memory allocation error\n");
	return(1);
      }
    }

  } else {

/*------- non template ----------*/

    ch = pfi->name;

  }

/*---
    Open this GRIB file and position to start of first record 
---*/


  if(diag) printf("22222 After open test: update= %d update_mape= %d\n",update,update_map);

  if(update_map==0) {

    if(!quiet)  printf("gribmap:  Opening GRIB file %s\n",ch);
    gfile = fopen(ch,"rb");

    if (gfile==NULL) {
      if (pfi->tmplat) {
	printf ("Warning: Could not open GRIB file: %s\n",ch);
	continue;
      } else {
	printf ("Could not open GRIB file.  File name is:\n");
	printf ("  %s\n",ch);
	return(1);
      }
    }

    if (pfi->tmplat) free(ch);
  
    /* Get file size */
  
    fseek (gfile,0L,2);
    flen = ftell(gfile);

    /*---
      Set up to skip appropriate amount and position 
      ---*/
  
    if (skip > -1) fpos = skip;
    else {
      fseek (gfile,0,0);
      rc = fread (rec,1,100,gfile);
      if (rc<100) {
	printf ("I/O Error reading header\n");
	return(1);
      }
      len = gagby(rec,88,4);
      fpos = len*2 + 100;
    }


    /*--- main loop 
      
      1) read a grib record (gribhdr)
      2) compare to each 2-d variable in the 5-D data volume
      defined by the dd for a match (gribrec)
      
      ---- main loop */

    irec=1;
    while (1) {
      rc=gribhdr(&ghdr);
      if (rc) break;
      rcgr=gribrec(&ghdr,pfi,pindx);
      if(rcgr==0) didmatch=1;
      if(rcgr>=100) didmatch=rcgr;
      irec++;
    }

    /*--- 
      see how we did
      ---*/

    if (rc==50) {
      printf ("I/O Error reading GRIB file\n");
      printf ("Possible cause: premature EOF\n");
      break;
    }
    if (rc>1 && rc!=98) {
      printf ("GRIB file format error \n");
      printf ("rc = %i\n",rc);
      return(rc);
    }

/*--- 
  break out if not templating
---*/

    if (!pfi->tmplat) break;
 
/*---
 end for update option
---*/

  } 

/*---
  return for another file if templating
---*/

}
/*--- 
  ALL DONE the last(only) GRIB file
---*/

if(!quiet) printf ("gribmap:  Reached EOF\n");

/*---
  check if file closed already for case where template was set,
  but it was not templated and the template code above closed it.
---*/

if(!gfile) fclose (gfile);

/*--- 
  open the map file 
---*/

if(mfile!=NULL) {
  if(diag) printf("-------------- file already open.....\n");
} else {
  if(write_map) {
    mfile = fopen(pfi->mnam,"wb");
  } else {
    printf("\nWARNING!!!! NOT writing the map !!!\n\n");
  }
}

if (mfile==NULL && write_map) {
  printf ("Could not open index file: %s\n",pfi->mnam);
  return(1);
} else if(write_map) {
  if(!quiet) printf("gribmap:  Writing the map...\n\n");
}

/*--- 
  output the map depending on version #
---*/

if(ver==1) {
  fwrite (pindx,sizeof(struct gaindx),1,mfile);
  if (pindx->hinum>0) fwrite(pindx->hipnt,sizeof(int),pindx->hinum,mfile);
  if (pindx->hfnum>0) fwrite(pindx->hfpnt,sizeof(float),pindx->hfnum,mfile);
  if (pindx->intnum>0) fwrite(pindx->intpnt,sizeof(int),pindx->intnum,mfile);
  if (pindx->fltnum>0) fwrite(pindx->fltpnt,sizeof(float),pindx->fltnum,mfile);
  fclose (mfile);
}

/*---
  new machine independent version
---*/

if(ver==2) {

  rc=wtgmap();
  if(rc == 601) {
    printf("\nGRIBMAP failed because of overflow in float -> IBM float conversion!!!\n");
    printf("Contact developer fiorino@llnl.gov...\n\n");
    return(601); 
  }
  fclose (mfile);
}

return(didmatch);

/*--- errors ---*/

 err8:
printf("Open Error:  Memory allocation Error\n");
return(99);

}

/* Read a GRIB header, and get needed info out of it.  */

int gribhdr (struct grhdr *ghdr) {
struct dt atim;
char rec[500],*ch,*gds;
int i,len ,rc,sign,mant;
int cpos;

  if (fpos+50>=flen) return(1);


    /* look for data between records BY DEFAULT */ 

    i = 0;
    fpos += i;
    rc = fseek(gfile,fpos,0);
    if (rc) return(50);
    ch=&rec[0];
    rc = fread(ch,sizeof(char),4,gfile);
    
    while ( (fpos < flen-4) && (i < scanlim) && !( *ch=='G' && *(ch+1)=='R' &&
			   *(ch+2)=='I' && *(ch+3)=='B' ) ) {
      fpos++;
      i++;
      rc = fseek(gfile,fpos,0);
      if (rc) return(50);
      rc = fread(ch,sizeof(char),4,gfile);
      if (rc<4) return(50);
    } 
    
    if (i == scanlim ) {
      printf ("GRIB header not found in scanning between records\n");
      printf ("->%c%c%c%c<-\n",*rec,*(rec+1),*(rec+2),*(rec+3));
      if(scaneof==1 || scanEOF==1) {
	if(scaneof) {
	  printf("...... punching out -e option set ......\n"); 
	  return(98);
	}
	if(scanEOF) {
	  printf("...... junk found, continue -E option set ......\n"); 
	  return(0);
	}
      } else {
	return(52);
      }
    } else if (fpos == flen -4) {
      if(scaneof==1 || scanEOF==1) {
	if(scaneof) {
	  printf("...... punching out -e option set flen-4 check ......\n"); 
	  return(98);
	}
	if(scanEOF) {
	  printf("...... junk found, continue -E option set flen-4 check ......\n"); 
	  return(0);
	}
      } else {
	return (53);
      }
    }

   /* SUCCESS redo the initial read */    
    
    rc = fseek(gfile,fpos,0);
    if (rc) return(50);
    rc = fread(rec,1,8,gfile);
    if (rc<8){
      if(fpos+8 >= flen) return(61);
      else return(62);
    }



  cpos = fpos;
  ghdr->vers = gagby(rec,7,1);
  if (ghdr->vers>1) {
    printf ("File is not GRIB version 0 or 1, 0 or 1 is required. \n");
    printf ("Version number is %i\n",ghdr->vers);
    if(scaneof) {
      printf("...... punching out -e option set junk at end ......\n"); 
      return(98);
    } else {
      return (99);
    }
  }

  if (ghdr->vers==0) {
    cpos += 4;
    rc = fseek(gfile,cpos,0);
    if (rc) return(50);
  } else {
    ghdr->len = gagby(rec,4,3);
    cpos = cpos + 8;
    rc = fseek(gfile,cpos,0);
    if (rc) return(50);
  }

  /* Get PDS length, read rest of PDS */

  rc = fread(rec,1,3,gfile);
  if (rc<3) return(50);
  len = gagby(rec,0,3);
  ghdr->pdslen = len;
  cpos = cpos + len;
  rc = fread(rec+3,1,len-3,gfile);
  if (rc<len-3) return(50);

  /* Get info from PDS */

  ghdr->id = gagby(rec,6,1);
  ghdr->gdsflg = gagbb(rec+7,0,1);
  ghdr->bmsflg = gagbb(rec+7,1,1);
  ghdr->parm = gagby(rec,8,1);
  ghdr->ltyp = gagby(rec,9,1);
  ghdr->level = gagby(rec,10,2);
  ghdr->l1 = gagby(rec,10,1);
  ghdr->l2 = gagby(rec,11,1);
  if (mpiflg) {        /* Use initial time from the descriptor for mpiflg */
    ghdr->btim.yr = *(pfi->abvals[3]);
    ghdr->btim.mo = *(pfi->abvals[3]+1);
    ghdr->btim.dy = *(pfi->abvals[3]+2);
    ghdr->btim.hr = *(pfi->abvals[3]+3);
    ghdr->btim.mn = *(pfi->abvals[3]+4);
  } else {
    ghdr->btim.yr = gagby(rec,12,1);
    ghdr->btim.mo = gagby(rec,13,1);
    ghdr->btim.dy = gagby(rec,14,1);
    ghdr->btim.hr = gagby(rec,15,1);
    ghdr->btim.mn = gagby(rec,16,1);
  }
  if (ghdr->btim.hr>23) ghdr->btim.hr = 0;  /* Special for NCAR */
  if (len>24) {
    ghdr->cent = gagby(rec,24,1);
    ghdr->btim.yr = ghdr->btim.yr + (ghdr->cent-1)*100;
  } else {
    ghdr->cent = -999;
    if (!(mpiflg) || !(mfcmn.fullyear)) {
      if (ghdr->btim.yr>49) ghdr->btim.yr += 1900;
      if (ghdr->btim.yr<50) ghdr->btim.yr += 2000;
    }
  }
  ghdr->ftu = gagby(rec,17,1);
  ghdr->tri = gagby(rec,20,1);
  if (ghdr->tri==10) {
    ghdr->p1 = gagby(rec,18,2);
    ghdr->p2 = 0;
  } else {
    ghdr->p1 = gagby(rec,18,1);
    ghdr->p2 = gagby(rec,19,1);
  }

  ghdr->fcstt = ghdr->p1;
  if (ghdr->tri>1 && ghdr->tri<6) ghdr->fcstt=ghdr->p2;
  atim.yr=0; atim.mo=0; atim.dy=0; atim.hr=0; atim.mn=0;
  if (ghdr->ftu==0) atim.mn = ghdr->fcstt;
  else if (ghdr->ftu==1) atim.hr = ghdr->fcstt;
  else if (ghdr->ftu==2) atim.dy = ghdr->fcstt;
  else if (ghdr->ftu==3) atim.mo = ghdr->fcstt;
  else if (ghdr->ftu==4) atim.yr = ghdr->fcstt;
  else ghdr->fcstt = -999;
/* mf */
/*  if notau != 0 then FORCE the valid DTG to be the base DTG */ 
  if(notau) ghdr->fcstt = -999 ;
/* mf */

/*  add the forecast time to the time of this grib field */

  if (ghdr->fcstt>-900) {
    timadd(&(ghdr->btim),&atim);
    ghdr->dtim.yr = atim.yr;
    ghdr->dtim.mo = atim.mo;
    ghdr->dtim.dy = atim.dy;
    ghdr->dtim.hr = atim.hr;
    ghdr->dtim.mn = atim.mn;
  } else {
    ghdr->dtim.yr = ghdr->btim.yr;
    ghdr->dtim.mo = ghdr->btim.mo;
    ghdr->dtim.dy = ghdr->btim.dy;
    ghdr->dtim.hr = ghdr->btim.hr;
    ghdr->dtim.mn = ghdr->btim.mn;
  }
  if (len>25) {
    ghdr->dsf = (float)gagbb(rec+26,1,15);
    i = gagbb(rec+26,0,1);
    if (i) ghdr->dsf = -1.0*ghdr->dsf;
    ghdr->dsf = pow(10.0,ghdr->dsf);
  } else ghdr->dsf = 1.0;


  /* If it is there, get info from GDS */

  if (ghdr->gdsflg) {
    rc = fread(rec,1,3,gfile);
    if (rc<3) return(50);
    len = gagby(rec,0,3);
    ghdr->gdslen = len;
    cpos = cpos + len;

/*mf --- malloc for my generic grid where I code the lon/lats in the GDS --mf*/

    gds=(char *)malloc(len+3);
    if(gds == NULL) {
      return(51);
    }

    rc = fread(gds+3,1,len-3,gfile);
    if (rc<len-3) return(50);
    ghdr->gtyp = gagby(gds,4,1);
    ghdr->gicnt = gagby(gds,6,2);
    ghdr->gjcnt = gagby(gds,8,2);
    ghdr->gsf1 = gagbb(gds+27,0,1);
    ghdr->gsf2 = gagbb(gds+27,1,1);
    ghdr->gsf3 = gagbb(gds+27,2,1);

/*
   printf ("gds len,typ,i,j,flags = %i %i %i %i %i %i %i\n",
   len,ghdr->gtyp,ghdr->gicnt,ghdr->gjcnt,ghdr->gsf1,ghdr->gsf2,ghdr->gsf3);
*/
    free(gds);

  } else ghdr->gdslen = 0;

  /* Get necessary info about BMS if it is there */

  if (ghdr->bmsflg) {
    rc = fread(rec,1,6,gfile);
    if (rc<6) return(50);
    len = gagby(rec,0,3);
    ghdr->bmsflg = len;
    ghdr->bnumr = gagby(rec,4,2);
    ghdr->bpos = cpos+6;
    cpos = cpos + len;
    rc = fseek(gfile,cpos,0);
  } else ghdr->bmslen = 0;

  /* Get necessary info from data header */

  rc = fread(rec,1,11,gfile);
  if (rc<11) return(50);
  len = gagby(rec,0,3);
  ghdr->bdslen = len;
  ghdr->iflg = gagbb(rec+3,0,2);
  i = gagby(rec,4,2);
  if (i>32767) i = 32768-i;
  ghdr->bsf = pow(2.0,(float)i);

  i = gagby(rec,6,1);
  sign = 0;
  if (i>127) {
    sign = 1;
    i = i - 128;
  }
  mant = gagby(rec,7,3);
  if (sign) mant = -mant;
  ghdr->ref = (float)mant * pow(16.0,(float)(i-70));

  ghdr->bnum = gagby(rec,10,1);
  ghdr->dpos = cpos+11;
/*
  printf ("len,flg,bsf,ref,bnum %i %i %g %g %i\n",len,ghdr->iflg,
    ghdr->bsf,ghdr->ref,ghdr->bnum);
*/

  if (ghdr->vers==0) {
    fpos = fpos + 8 + ghdr->pdslen + ghdr->gdslen +
                      ghdr->bmslen + ghdr->bdslen;
  } else fpos = fpos + ghdr->len;

  return(0);

}

/* Routine to determine the location of the GRIB record in terms of
   the GrADS data set, and fill in the proper values at the proper
   slot location. */

int gribrec (struct grhdr *ghdr, struct gafile *pfi, struct gaindx *pindx) {
float (*conv) (float *, float);
struct gavar *pvar;
int i,ioff,iz,it,joff,nsiz,flag;
float z,t;

  /* Verify that we are looking at the proper grid type */

  joff =0;
  nsiz = nrec * nelements ;
  if (ghdr->iflg) {
    if (verb) {
      printf ("GRIB record contains harmonic or complex packing\n");
      printf ("  Record is skipped.\n");
      printf ("  Variable is %i\n",ghdr->parm);
    }
    return(10);
  }
  if (pfi->grbgrd==255 || pfi->grbgrd<-900) {
    if (!ghdr->gdsflg) {
      if (verb) {
        printf ("GRIB record contains pre-defined grid type: "); 
        printf ("GrADS descriptor specifies type 255\n");
	gribpr(ghdr);
      }
      return(20);
    } 
    if ( pfi->ppflag) {
      if ( (ghdr->gicnt != pfi->ppisiz) || (ghdr->gjcnt != pfi->ppjsiz) ) {
        if (verb) {
          printf ("GRIB grid size does not match descriptor: "); 
	  gribpr(ghdr);
        }
        return(300);
      }
    } else {
      if ( (ghdr->gicnt != pfi->dnum[0]) || (ghdr->gjcnt != pfi->dnum[1]) ) {
        if (verb) {
          printf ("GRIB grid size does not match descriptor:");
	  gribpr(ghdr);
        }
        return(301);
      }
    }
  } else {

/*--- 
  special case for 
  GRIB grid number (dtype grib NNN) == 29
---*/

    if (pfi->grbgrd==29) {
      if (ghdr->id!=29 && ghdr->id!=30) {
        if (verb) {
	  printf("Record has wrong GRIB grid type: ") ; 
	  gribpr(ghdr);
	}
        return(400);     
      }
      if (ghdr->id==29) joff = nelements;
      nsiz = 2 * nelements ;
    } else {
      if (ghdr->id != pfi->grbgrd) {
        if (verb) {
	  printf("%s","Record has wrong GRIB grid type: "); 
	  gribpr(ghdr);
	}
        return(401);     
      }
    }
  }
    
  /* Calculate the grid time for this record.  If it is non-integer
     or if it is out of bounds, just return. */

/* only look for tauoff hour taus */

  if( tauflg && (ghdr->ftu==1 && ghdr->fcstt!=tauoff) ) {
    if (verb) {
      printf("%s %d","--f-- Forecast Time does not match : ",tauoff);
      gribpr(ghdr);
    }
    return(32);
  }

/*
  printf("qqq %i%i%i%i%i ... %i%i%i%i%i\n",
       ghdr->btim.yr,ghdr->btim.mo,ghdr->btim.dy,ghdr->btim.hr,ghdr->btim.hr,
       btimdd.yr,btimdd.mo,btimdd.dy,btimdd.hr,btimdd.mn);
*/

/*  check if base DTG in grib field matches initial time in DD */

  if(tau0 &&
     ((ghdr->btim.yr != btimdd.yr) ||
     (ghdr->btim.mo != btimdd.mo) ||
     (ghdr->btim.dy != btimdd.dy) ||
     (ghdr->btim.hr != btimdd.hr) ||
     (ghdr->btim.mn != btimdd.mn)) ) {
    if (verb) {
      printf("%s","--b-- Base Time does not match Initial Time in DD: "); 
      gribpr(ghdr);
    }
    return(34);
  }

  t = t2gr(pfi->abvals[3],&(ghdr->dtim));

  if (t<0.99 || t>((float)(pfi->dnum[3])+0.01)) {
    if (verb) {
      printf("%s","----- Time out of bounds: "); 
      gribpr(ghdr);
    }
    return(36);
  }

  it = (int)(t+0.01);

  if (fabs((float)it - t)>0.01) {
    if (verb) {
      printf("----- Time non-integral. %g %g: ",(float)it,t);  
      gribpr(ghdr);
    }
    return(38);
  }
/*  it = (int)(t+0.1); */
  it = (it-1)*pfi->trecs;

  /* See if we can match up this grid with a variable in
     the data descriptor file */

  pvar = pfi->pvar1;
  i = 0;
  flag=0;


while (i<pfi->vnum) {

  if (pvar->levels>0) {

/*
  multi level data
*/

    if (pvar->units[0]==ghdr->parm && pvar->units[1]==ghdr->ltyp) {
    /*
      mf 940602 - look for time range indicator match
    */
      if(diag) printf("qqq1 %d %d %d %d\n",pvar->units[3],ghdr->tri,pvar->units[0],pvar->units[1]); 
      if(pvar->units[3] < -900 || pvar->units[3] == ghdr->tri) {
        conv = pfi->ab2gr[2];
        z = conv(pfi->abvals[2],ghdr->level);
	if(diag) printf("qqq2 %g %g \n",z,(float)(pvar->levels)); 
        if (z>0.99 && z<((float)(pvar->levels)+0.01)) {
          iz = (int)(z+0.5);
	  /* 
	    search for level match
	    */
	  if(diag) printf("qqq3 %d\n",iz); 
          if ( fabs(z-(float)iz)<0.01) {
            iz = (int)(z+0.5);
	    if(diag) printf("qqq4 match\n"); 
            ioff = pvar->recoff + iz - 1;
            gribfill (it+ioff,joff,nsiz,ghdr,pindx);
            flag=1;
            i = pfi->vnum + 1;   /* Force loop to stop */
          }
        }
      }
    }
  } else {

/* 
  sfc data
*/
    if (pvar->units[0]==ghdr->parm && pvar->units[1]==ghdr->ltyp) {
      if ( (pvar->units[3]<-900 && pvar->units[2]==ghdr->level) ||
	  (pvar->units[3]>-900 && pvar->units[2]==ghdr->l1 && pvar->units[3]==ghdr->l2 )

/*mf 940602
  test for time range indicator
  NB!!!!  - if a variable has all the same parameters except the
  time range indicator then explicitly set it up!!
*/

/*mf 951129	

  took this test out.... require match in tri test for level (2 byte number) only
  || (pvar->units[3] == ghdr->tri && pvar->units[2]==ghdr->l1)

*/
	  || (pvar->units[3] == ghdr->tri && pvar->units[2]==ghdr->level) ) {
	ioff = pvar->recoff;
	gribfill (it+ioff,joff,nsiz,ghdr,pindx);
	i = pfi->vnum+1;  /* Force loop to stop */
	flag=1;
      }
    }
  }
  pvar++; i++;
}


if (flag && verb) printf("!!!!! MATCH: "); 
if (!flag && verb) printf("..... NOOOO: "); 
if (verb) gribpr(ghdr); 

if(flag) return(0);
if(!flag) return(1);

}


/* Fill in values for this record, now that we have found how it
   matches.  We are not handling the time aspect as yet */

void gribfill (int ioff, int joff, int nsiz, struct grhdr *ghdr,
      struct gaindx *pindx) {
  ioff = nsiz*ioff + joff;
  *(pindx->intpnt+ioff) = ghdr->dpos;
  if (ghdr->bmsflg) *(pindx->intpnt+ioff+1) = ghdr->bpos;
  *(pindx->intpnt+ioff+2) = ghdr->bnum;
  *(pindx->fltpnt+ioff) = ghdr->dsf;
  *(pindx->fltpnt+ioff+1) = ghdr->bsf;
  *(pindx->fltpnt+ioff+2) = ghdr->ref;
}

void gribpr(struct grhdr *ghdr) {

  printf ("% 4i % 10i % 3i % 1i % 1i % 3i % 3i %-5i % 10i % 10i % 2i ",
    irec,fpos,ghdr->id,ghdr->gdsflg,ghdr->bmsflg,ghdr->parm,ghdr->ltyp,
    ghdr->level,ghdr->dpos,ghdr->bpos,ghdr->bnum);
  printf ("btim: %04i%02i%02i%02i:%02i ",ghdr->btim.yr,
    ghdr->btim.mo,ghdr->btim.dy,ghdr->btim.hr,ghdr->btim.mn);
  printf ("tau: % 6i ",ghdr->fcstt);
  printf ("dtim: %04i%02i%02i%02i:%02i ",ghdr->dtim.yr,
    ghdr->dtim.mo,ghdr->dtim.dy,ghdr->dtim.hr,ghdr->dtim.mn);
  printf("\n");

}

int wtgmap(void) {

int i,nb,bcnt,idum;
float fdum;
unsigned char *map;
unsigned char ibmfloat[4];

/* calculate the size of the ver==1 index file */

nb=
  2 + (4*4) +                          /* version in byte 2, then 4 ints with number of each data type */
    pindx->hinum*sizeof(int)+
      pindx->hfnum*sizeof(int)+
	pindx->intnum*sizeof(int)+
	  pindx->fltnum*sizeof(float) ;

/* add additional info */

nb=nb+7;      /* base time (+ sec)  for compatibility with earlier version 2 maps */
nb=nb+8*4;    /* grvals for time <-> grid conversion */

if(diag) {
  printf("\nqqq nb= %d\n",nb);
  printf("    btimdd.yr = %d\n",btimdd.yr);
  printf("    btimdd.mo = %d\n",btimdd.mo);
  printf("    btimdd.dy = %d\n",btimdd.dy);
  printf("    btimdd.hr = %d\n",btimdd.hr);
  printf("    btimdd.mn = %d\n",btimdd.mn);
}
/* allocate space for the map */

map = (unsigned char *)malloc(nb);
if(map == NULL) {
  fprintf(stderr,"malloc error creating the map in gagmap\n");
  return(60);
}


if(diag) {

  for(i=0;i<pindx->hinum;i++) {
    printf("i = %d pindx->hipnt = %d\n",i,*(pindx->hipnt+i));
  }
  for(i=0;i<pindx->hfnum;i++) {
    printf("i = %d pindx->hfpnt = %g\n",i,*(pindx->hfpnt+i));
  }
  for(i=0;i<pindx->intnum;i++) {
    printf("i = %d pindx->intpnt = %d\n",i,*(pindx->intpnt+i));
  }
  for(i=0;i<pindx->fltnum;i++) {
    printf("i = %d pindx->fltpnt = %g\n",i,*(pindx->fltpnt+i));
  }

}
bcnt=0;

gapby(0,map,bcnt,1);             bcnt++  ;   /* set the first byte to 0 */
gapby(ver,map,bcnt,1);           bcnt++  ;   /* set the second byte to the version number */

putint(pindx->hinum,map,&bcnt);              /* # ints in header   */
putint(pindx->hfnum,map,&bcnt);              /* # floats in header   */
putint(pindx->intnum,map,&bcnt);             /* # index ints   */
putint(pindx->fltnum,map,&bcnt);             /* # index floats   */

/*--- 
   write out base time for consistency with earlier ver 2 maps
   not really needed
---*/

gapby(btimdd.yr,map,bcnt,2);  bcnt+=2 ;   /* initial year */
gapby(btimdd.mo,map,bcnt,1);  bcnt++  ;   /* initial month */ 
gapby(btimdd.dy,map,bcnt,1);  bcnt++  ;   /* initial day */
gapby(btimdd.hr,map,bcnt,1);  bcnt++  ;   /* initial hour */
gapby(btimdd.mn,map,bcnt,1);  bcnt++  ;   /* initial minute */
gapby(0,map,bcnt,1);          bcnt++  ;   /* initial second */


/*---
  load header
---*/

if(pindx->hinum) {

  for(i=0;i<pindx->hinum;i++) {
    idum=*(pindx->hipnt+i);
    putint(idum,map,&bcnt);
    if(diag) printf("i = %d pindx->hipnt = %d\n",i,idum);
  }

}

if(pindx->hfnum) {

/* blank for now */

}

/*---
  load index
---*/

for(i=0;i<pindx->intnum;i++) {
  idum=*(pindx->intpnt+i);
  putint(idum,map,&bcnt);
  if(diag) printf("i = %d pindx->intpnt = %d\n",i,idum);
}

for(i=0;i<pindx->fltnum;i++) {
  fdum=*(pindx->fltpnt+i);
  rc=flt2ibm(fdum, ibmfloat); 
  if(rc<0) return(601);
  if(diag) printf("i = %d pindx->fltpnt = %g\n",i,fdum);
  strncpy(&map[bcnt],ibmfloat,4); bcnt+=4;
}


/* write out the factors for converting from grid to absolute time */ 

for(i=0;i<8;i++) {
  fdum=*(pfi->grvals[3]+i);
  rc=flt2ibm(fdum, ibmfloat); 
  if(rc<0) return(601);
  if(diag) printf("i = %d grvals = %g\n",i,fdum);
  strncpy(&map[bcnt],ibmfloat,4); bcnt+=4;
}


if(diag) printf("qqqqq the size of the map is %d the calc size is %d\n",bcnt,nb);

if(write_map) {
  fwrite(map,1,bcnt,mfile);
}
free(map);

return(0);

}
/*---

  routine to dump a 4 byte int into a character stream 
  Mike Fiorino 960829

---*/

void putint(int dum, unsigned char *buf,int *off) {

  int offset;
  offset=*off;
  if(dum < 0) {
    dum=-dum;
    gapby(dum,buf,offset,4);
    gapbb(1,buf,offset*8,1);
  } else {
    gapby(dum,buf,offset,4);
  }
  offset+=4;
  *off=offset;

}
