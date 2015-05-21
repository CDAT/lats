/*
   newdtg:  Return a new dtg value by adding/subtracting the a given
            dtg value with some increment value.

   options: 

*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>

struct cxclock {
  time_t year ;
  time_t month ;
  time_t date ;
  time_t hour ;
  time_t minute ;
  time_t second ; 
  char timezone[32] ;
  char clockenv[64] ;
  int julian_day ;
  time_t epoch_time_in_sec ;
} ;

/* global variable of the time object */
static struct cxclock timeobj;

void get_tm_struct(time_t t)
{
  struct tm *ptr_tm ;
  char buf[10] ;
  time_t sec ;
  sec=t ;

  /* get tm struct from seconds */
#ifdef SUNOS
  ptr_tm=gmtime(&sec) ;
#else
  ptr_tm=localtime(&sec) ;
#endif

  if ( ! ptr_tm ) {
    fprintf(stderr,
    "error in gmtime function,errno = %d\n",errno) ;
    exit(-1) ;
  }
	
  /* copy tm struct to cxclock private members */

  timeobj.epoch_time_in_sec=sec ;
  timeobj.year=ptr_tm->tm_year + 1900  ;
  timeobj.month=(ptr_tm->tm_mon) + 1 ;
  timeobj.date=ptr_tm->tm_mday ;
  timeobj.hour=ptr_tm->tm_hour ;
  timeobj.minute=ptr_tm->tm_min ;
  timeobj.second=ptr_tm->tm_sec ;

#ifdef SUNOS
  strncpy(timeobj.timezone,
          ptr_tm->tm_zone,sizeof(timeobj.timezone)-1) ;   
#else
  /* 
     tzname is an external variable which contain 
     timezone string for SVR5 Unix system 
     (SGI, Solaris,HP,...) 
  */
  strncpy(timeobj.timezone,
          tzname[0], sizeof(timeobj.timezone)-1) ;
#endif

  /* get julian date */
  strftime(buf, sizeof(buf),"%j", ptr_tm) ;
  timeobj.julian_day=atoi(buf) ;
}

void sys_time()
{
  time_t sec ;
  time(&sec) ; /* get the GMT offset from the local machine */
  if (sec == -1 ) {
    fprintf(stderr,
            "error in time function,errno = %d\n"
            ,errno) ;
    exit(-1) ;
  }

  /* Set tm struct using new epoch time */
  get_tm_struct(sec) ;
  return;
}

void advance(const int offset_hrs)
{
#define SEC_TO_HR 3600
  /* Compute new seconds form epoch */
  int offset_sec ;
  time_t new_epoch_time	;
  offset_sec=offset_hrs * SEC_TO_HR ;
  /* printf("offset_sec %d %d \n",offset_sec,
          timeobj.epoch_time_in_sec ); */

  new_epoch_time = timeobj.epoch_time_in_sec + 
                   offset_sec ;   

  /* printf("new_epoch_time =  %d \n",
        new_epoch_time ); */
	
  /* Set time using new epoch time */
  get_tm_struct(new_epoch_time) ;
  return;
}

void dtg_parser(const char* dtg_input)
{
  char tmpstr[5] ;

  if ( strlen(dtg_input) != 10 )
  {
    fprintf(stderr,"Day Time Group format error: \n");
    fprintf(stderr,"\tThe correct format is YYYYMMDDHH\n"); 
    fprintf(stderr,"\twhile your input DTG string is %s\n", 
            dtg_input); 
    fprintf(stderr,"\t......aborting\n") ;
    exit(-1) ;
  }

	tmpstr[0]=dtg_input[0] ;
	tmpstr[1]=dtg_input[1] ;
	tmpstr[2]=dtg_input[2] ;
	tmpstr[3]=dtg_input[3] ;
	tmpstr[4]='\0' ;
	timeobj.year=atoi(tmpstr) ;

	tmpstr[0]=dtg_input[4] ;
	tmpstr[1]=dtg_input[5] ;
	tmpstr[2]='\0' ;
	timeobj.month=atoi(tmpstr) ;

	tmpstr[0]=dtg_input[6] ;
	tmpstr[1]=dtg_input[7] ;
	tmpstr[2]='\0' ;
	timeobj.date=atoi(tmpstr) ;

	tmpstr[0]=dtg_input[8] ;
	tmpstr[1]=dtg_input[9] ;
	tmpstr[3]='\0';   
	timeobj.hour=atoi(tmpstr) ;
	
}

void set_user_time(const char* s)
{
  time_t sec ;
  char buf[10] ;
  struct tm *ptr_tm ;
  dtg_parser(s) ;

  time(&sec) ; /* get the GMT offset from the local machine */
  if (sec == -1 ) {
    fprintf(stderr,
    "error in time function,errno = %d\n",errno) ;
    exit(-1) ;
  }
	
  ptr_tm=gmtime(&sec) ;
  if ( ! ptr_tm ) {
    fprintf(stderr,
            "error in gmtime function,errno= %d\n"
            ,errno) ;
    exit(-1) ;
  }

  /* convert the DTG date to seconds */
  ptr_tm->tm_hour=(int) timeobj.hour ; 
  ptr_tm->tm_mday=(int) timeobj.date  ; 
  ptr_tm->tm_mon=(int) (timeobj.month-1) ; 
  ptr_tm->tm_year=(int)(timeobj.year -1900 );  

  /* get julian date
     Convert a tm structure to a calendar time 
     (seconds since 00:00:00 UTC, January 1, 1970).
     */

#ifdef SUNOS
   sec=timegm(ptr_tm) ;
#else
   sec=mktime(ptr_tm) ;
#endif

   ptr_tm=gmtime(&sec) ;
   timeobj.epoch_time_in_sec=sec ;

   strftime(buf, sizeof(buf),"%j", ptr_tm) ;
   timeobj.julian_day=atoi(buf) ;
   timeobj.minute=0 ;
   timeobj.second=0 ;
   return;
}

void reset(const char* s)
{
  set_user_time(s) ;
  return;
}

/* Function to convert DTG string to 
   GrADS Time format
   */

char* dtg2grads_time(char *dtg) ;
int is_signed_integer(char *str) ;

int main(int argc , char* argv[])
{
   static int offset = 0  ;
   static char tzone[]="TZ=GMT" ;
   char dtgstr[32] ;

   int i ;
   /* Process command line arguments */
   i=1 ;
   if ( argc != 3 ) {
     fprintf(stdout,"Usage: %s current-dtg (YYYYMMDDHH) increment(HH) \n",
             argv[0]) ;
     return 1 ;
   } else {
      if ( strlen(argv[1]) != 10 ) {
         fprintf(stderr,"DTG must be in YYYYMMDDHH format.\n") ;
         return 1 ;
      } else {
      	strcpy(dtgstr,argv[1]) ;
      }
      /* Any valid integer value is taken as 
       an offset to current DTG value */
      if ( is_signed_integer(argv[2]) == 1 )
      {    
         offset=atoi(argv[2])  ;
         /* printf("offset= %d\n",offset); */
      } else {
         fprintf(stderr,
                 "Increment value: %s is not a valid integer.\n",argv[2]) ;
         return 1 ;
      }
   }

   /* 
      For Operational DTG, we always use GMT time
      Set TZ to GMT before instantiate cxclock object
   */
 
   putenv(tzone) ;

   sys_time() ;    

   /*
     Set the gmt time to the current clock value to
     CRDATE or current DTG TIME.
     */
   reset(dtgstr) ;   

   /* Adjust the value if we have any offset */
   if ( offset != 0 ) 
   {
     /*
       Adjust the clock value by the offset value
       */
      advance(offset) ;
      sprintf(dtgstr,"%04d%02d%02d%02d\0",
              timeobj.year,timeobj.month,
	      timeobj.date,timeobj.hour) ;
   }
   fprintf(stdout,"%s\n",dtgstr) ;
    
return 0 ;

}


/* This function return 0 or 1 depends on if the 
   str is a valid integer value */
int is_signed_integer(char *str)
{
   int nlen, i  ; 
   int beg ;
   if ( str==NULL ) return 0 ;

   nlen=strlen(str) ;

   /* Check the first character of str 
      (look for - or + sign)
      */
   if ( str[0] == '-' || str[0] == '+' ) 
       beg = 1 ;
   else 
       beg = 0 ;

   for ( i = beg ; i < nlen ; i++ )
   {
      if ( ! isdigit(str[i]) ) return 0 ;
   }

   return 1 ;
}


