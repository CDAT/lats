#define _POSIX_SOURCE 1
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lats.h"
#include "latsint.h"
#include "latsparm.h"

char *_lats_routine_name_ = (char *)0;
void command_line_help(char *) ;

main (int argc, char *argv[])  {

int i,j;
int help=0;
int vertdim=0;
int center=0;
int parm=0;
int delim=0;
int all=0;
int table=0;
int header=1; 

int cmdflg = 0;
int hstflg = 1; /*mf default is to USE command line editing mf*/
int gflag = 0;
char *icmd = NULL;

char rec[256],rec1[10],rec2[10];
char fmt[256];

if (argc>1) {
  for (i=1; i<argc; i++) {
    if (*(argv[i])=='-' && *(argv[i]+1)=='h' && *(argv[i]+2)=='e' && *(argv[i]+3)=='l' && *(argv[i]+4)=='p' ) {
      command_line_help(argv[0]);
      return(0);
    } else if (*(argv[i])=='-') {
      j = 1;
      while (*(argv[i]+j)) {
	if (*(argv[i]+j)=='v') vertdim = 1;
	else if (*(argv[i]+j)=='c') center = 1;
	else if (*(argv[i]+j)=='a') all = 1;
	else if (*(argv[i]+j)=='t') table = 1;
	else if (*(argv[i]+j)=='n') header=0;
	else if (*(argv[i]+j)=='p') parm = 1;
	else if (*(argv[i]+j)=='d') delim = 1;
	else printf ("Unknown command line option: %c\n",*(argv[i]+j));
	j++;
      }
    } else printf ("Unknown command line keyword: %s\n",argv[i]);
  }

} else {

  command_line_help(argv[0]);
  return(0);

} 

if(argc==2 && delim && j == 2 ) { 
  command_line_help(argv[0]);
  return(0);
}

if(!latsParmCreateDefaultTable()) return(1);

/*---
  set up
---*/

if(table || delim) header=0;

/*---
  parameter table
---*/

if(parm | all) {
 

  if(delim) {
    strcpy(fmt,"%s,%s,%s,%s,%d,%d,%d,%d\n");
  } else if(table) {
    strcpy(fmt,"%-8s | %3d | %-60s | %-10s | %-5s | %-7s | %4d | %4d | |\n");
  } else {
    strcpy(fmt,"%-8s | %-45s | %-10s | %-7s | %3d | %4d | %4d \n");
  }


  if(table) printf("#!variable\n");


  if(header) {
    printf("\n... LATS PARAMETERS ...\n\n");
    printf(
"name     | description                                   | units      | levtype | GRIB ID        \n");

    printf(
"                                                                                      | GRIB dec scale factor\n");
    printf(
"                                                                                             | GRIB fixed prec\n");
  }


  for(i=0;i<LATS_DEFAULT_NPARMS;i++) {

    if(delim || table) {
      strcpy(rec,latsDefaultParms[i].title);
    } else {
      strncpy(rec,latsDefaultParms[i].title,45);
    }

    if(table) {

      if(latsDefaultParms[i].datatype == LATS_FLOAT) strcpy(rec1,"float");
      if(latsDefaultParms[i].datatype == LATS_INT) strcpy(rec1,"int");

      printf(fmt,
	     latsDefaultParms[i].name,
	     latsDefaultParms[i].id,
	     rec,
	     latsDefaultParms[i].units,
	     rec1,
	     latsDefaultParms[i].levelname,
	     latsDefaultParms[i].scalefac,
	     latsDefaultParms[i].precision
	     );


    } else {

      printf(fmt,
	     latsDefaultParms[i].name,
	     rec,
	     latsDefaultParms[i].units,
	     latsDefaultParms[i].levelname,
	     latsDefaultParms[i].id,
	     latsDefaultParms[i].scalefac,
	     latsDefaultParms[i].precision
	     );

    }

  }

}

/*---
  vertical dimension table
---*/

if(vertdim | all) {

  if(delim) {
    strcpy(fmt,"%s,%s,%s,%s,%s,%d,%d,%d,%d\n");
  } else if(table) {
    strcpy(fmt,"%-8s | %-35s | %-10s | %s | %s | %d | %d | %d | %d \n");
  } else {
    strcpy(fmt,"%-8s | %-35s | %-10s | %-6s | %-4s | %3d | %3d | %4d | %5d \n");
  }

  if(header) {
    printf("\n... LATS VERTICAL DIMENSIONS ...\n\n");
    printf(
"name     | description                         | units      | verticality                             \n");
    printf(
"                                                                     | postive orientation            \n");
    printf(
"                                                                            | GRIB ID                 \n");
    printf(
"                                                                                  | GRIB PDS 11       \n");
    printf(
"                                                                                        | PDS12| PDS 11+12\n");
  }

  if(table) printf("#!vert\n");

  for(i=0;i<LATS_DEFAULT_NVERTS;i++) {

    if(delim) {
      strcpy(rec,latsDefaultVerts[i].descrip);
    } else {
      strncpy(rec,latsDefaultVerts[i].descrip,35);
    }
    

    if(latsDefaultVerts[i].verticality == LATS_SINGLE_LEVEL) strcpy(rec1,"single");
    if(latsDefaultVerts[i].verticality == LATS_MULTI_LEVEL) strcpy(rec1,"multi");
    if(latsDefaultVerts[i].positive == LATS_UP) strcpy(rec2,"up");
    if(latsDefaultVerts[i].positive == LATS_DOWN) strcpy(rec2,"down");

    printf(fmt,

	   latsDefaultVerts[i].name,
	   rec,
	   latsDefaultVerts[i].units,
	   rec1,
	   rec2,
	   latsDefaultVerts[i].gribid,
	   latsDefaultVerts[i].grib_p1,
	   latsDefaultVerts[i].grib_p2,
	   latsDefaultVerts[i].grib_p3
	     
	   );
    }
  
}

/*---
  center table
---*/

if(center | all) {

  if(table) printf("#!center\n");

  if(header) {
    printf("\n... LATS CENTERS ...\n\n");
    printf(
"name         | GRIB process ID\n");
    printf(
"                   | GRIB center ID\n");
    printf(
"                         | GRIB subcenter ID\n");
  }

  if(delim) {
    strcpy(fmt,"%s,%d,%d,%d\n");
  } else if(table) {
    strcpy(fmt,"%-12s | %3d | %3d | %3d \n");
  } else {
    strcpy(fmt,"%-12s | %3d | %3d | %3d \n");
  }

  for(i=0;i<LATS_DEFAULT_NCENTERS;i++) {

    printf(fmt,
	   latsDefaultCenters[i].center,
	   latsDefaultCenters[i].gribid,
	   latsDefaultCenters[i].grib_center,
	   latsDefaultCenters[i].grib_subcenter
	   );
  }
  
}


}

void command_line_help(char *name) {
/*--- 
  output command line options 
---*/

printf("%s -- Version " LATSLS_VERSION "\n\n",name);
printf("list the internal tables in the LATS library\n\n");
printf("Command line options: \n\n");
printf("          -help   Just this help\n");
printf("          -a      all tables\n");
printf("          -p      parameter table\n");
printf("          -c      center table\n");
printf("          -v      vertical dimension table\n");
printf("          -t      table format output\n");
printf("          -n      no descriptive headers\n");
printf("          -d      comma deliminated output\n\n");
printf("Options may be combined, e.g., -pcd\n\n");
printf("   Example:\n\n");
printf("   %s -ad \n\n",name);
printf("     list all internal LATS tables in comma deliminated format\n\n");

}
