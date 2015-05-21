      program testfh
C         
C         Sample program to test LATS_FORECAST_HOURLY
C         to demonstrate a typical AMIP II application
C         
C         Mike Fiorino, PCMDI (lats@pcmdi.llnl.gov)
C         version 1.0  
C         10 March, 1997
C         
C*******  WARNING 
C*******  machine / compiler dependent include line
C*******
C
C*******  there is no standard but this worked on all plaformts.
C*******  NB - this is a relative path, and needs to be adjusted
C*******  depending on where you installed LATS
C          
C         use .inc from development dir
C

      include "lats.inc"
C         
C         define the grid:
C         
C         72 points in long (no wrap)
C         46 points in lat (pole points)
C         3 levels
C         2 variables
C
      parameter(ni=72,nj=46,nk=3,nv=2)

      character*20 center
      character*20 model
      character*9 var

      dimension var(nv),id_var(nv)
      double precision rlon(ni),rlat(nj),plev(nk),slev
      real ta(ni,nj,nk),psl(ni,nj)
C         
C         define the production center and the process (model)
C 
      data center/'PCMDI'/
      data model/'lats'/
C         
C         set the variable names (AMIP II convention)
C         
C         psl - sea level pressure (Pa)
C         ta - air temperature (degK)
C
      data var/'psl','ta'/ 
C         
C         define the pressure levels (hPa)
C
      data plev/850.0,500.0,200.0/

C         
C         set forecast time parameters
C
C         idt - increment in hours between forecasts
C         iyr0,im0,ida0,ihr0 - base time
C         

      idt=3
      iyr0=1997	
      imo0=5
      ida0=17
      ihr0=0


C0000000000000000000000000000000000000000000000000
C         
C         STEP #0  --  set the parameter table
C
C0000000000000000000000000000000000000000000000000

cc      id_parmtab=latsparmtab("../table/fnmoc.lats.table")
cc      if(id_parmtab.eq.0) stop 'latsparmtab error'

C1111111111111111111111111111111111111111111111111
C         
C         STEP #1  -- create a LATS file defining 
C
C         1) convention (e.g., GRIB with a standard calandar)
C         2) time statistic (e.g., monthly)
C         3) time increment (1 for 1 month) 
C         4) who or the center (e.g., PCMDI)
C         5) what or the model (e.g., lats)
C
C1111111111111111111111111111111111111111111111111

      nconv=2
      do iconv=1,nconv

        if(iconv.eq.1) latsconv=LATS_GRADS_GRIB
        if(iconv.eq.2) latsconv=LATS_COARDS

        if(latsconv.eq.LATS_GRADS_GRIB) then
          id_fil = latscreate('fnmoc.fh',
     $         latsconv,
     $         LATS_STANDARD,
     $         LATS_FORECAST_HOURLY,idt,center,
     $         model,'LATS GRIB test')
          print*,'LATS grib file id = ',id_fil
        else if(latsconv.eq.LATS_COARDS) then
          id_fil = latscreate('fnmoc.fh',
     $         LATS_COARDS,
     $         LATS_STANDARD,
     $         LATS_FORECAST_HOURLY,idt,center,
     $         model,'LATS netcdf test')
          print*,'LATS netcdf file id = ',id_fil
          
        else
          go to 800
        endif



C1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a         
C
C         STEP #1a  -- set a basetime:
C
C         1) required for using LATS_FORECAST_HOURLY
C
C         do this only once !!!
C
C1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a1a         

        
        irc=latsbasetime(id_fil,iyr0,imo0,ida0,ihr0)
        print*,'irc from set latsbasetime = ',irc

C222222222222222222222222222222222222222222222222222222         
C
C         STEP #2  -- define the grid:
C
C         1) lons and lats and # points
C         2) type (e.g., linear, gaussian or generic)
C
C         do this only once !!!
C
C222222222222222222222222222222222222222222222222222222         

        if(iconv.eq.1) then

          do i=1,ni
            rlon(i)=0.0+(i-1)*360.0/ni
          end do

          do j=1,nj
            rlat(j)=-90.0+(j-1)*180.0/(nj-1)
          end do

          id_grd=latsgrid("u54",LATS_LINEAR, ni, rlon, nj, 
     $         rlat)
          if(id_grd.eq.0) stop 'latsgrid error'

        endif

C333333333333333333333333333333333333333333333333333333         
C
C         STEP #3  -- define the vertical dimension:
C
C         1) multi or sfc
C         2) levels
C
C333333333333333333333333333333333333333333333333333333         

        id_lev=latsvertdim('pressure', 'plev', nk, plev)
        if(id_lev.eq.0) stop 'latsvertdim error'

C444444444444444444444444444444444444444444444444444444         
C         
C         STEP #4  -- define the variables
C         
C         1) variable (from the table, e.g., "psl" for mean sea leve pressure)
C         2) data type (typically LATS_FLOAT)
C         3) time statistic type (e.g., LATS_AVERAGE)
C         4) associated vertical dimension from step #3
C
C
C444444444444444444444444444444444444444444444444444444         

        id_var(1)=latsvar(id_fil,var(1),
     $       LATS_FLOAT,LATS_INSTANT,id_grd,
     $       0, 'sfc variable')
        if(id_var(1).eq.0) stop 'latsvar(1) error'

        id_var(2)=latsvar(id_fil,var(2),
     $       LATS_FLOAT,LATS_INSTANT,id_grd,
     $       id_lev, 'ua variable')
        if(id_var(2).eq.0) stop 'latsvar(2) error'
	
C4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a 
C
C         STEP #4a (optional) - set a missing value 
C         uncomment the following line
C
CCCC      ierr=latsmissreal(id_fil,id_var(1),1e20,1e13)
C
C4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a4a 

C555555555555555555555555555555555555555555555555555555         
C         
C         STEP #5 -- get and write the fields
C
C         1) define the valid time (beginning of the interval)
C         2) read(generate) some data
C         3) write it out
C
C555555555555555555555555555555555555555555555555555555         

C         
C         valid time 00UTC 1 Jan 1979
C         
        iyr=0
        imo=0
        ida=0

        do ihr=0,6,idt
C         
C         create sfc field and write
C
	
          call read_data(psl,var(1),ni,nj,0,imo)
C
C	I know this is overkill, but if you don't do this the hp version
C	will core dump!
C
	  slev=0.0D0
          id_write1=latswrite(id_fil,id_var(1),slev,
     $         iyr,imo,ida,ihr,psl)
          if(id_write1.eq.0) stop 'latswrite error - sfc'
C         
C         create multi-level field
C
          do k=1,nk
            call read_data(ta(1,1,k),var(2),ni,nj,k,imo)
            id_write2=latswrite(id_fil,id_var(2),plev(k),
     $           iyr,imo,ida,ihr,
     $           ta(1,1,k))
            if(id_write2.eq.0) stop 'latswrite error - ua'
          end do
          
        end do

C666666666666666666666666666666666666666666666666666666         
C         
C         STEP #6 - close the file
C         
C         1) for the GRIB, this creates the GrADS .ctl file for
C         futher processing by cdunif, GrADS or VCS
C
C666666666666666666666666666666666666666666666666666666         

        id_close = latsclose(id_fil)
        if(id_close.eq.0) stop 'latsclose error'

      end do
C         
C         normal exit
C         
      go to 999

C         
C         error conditions
C         
 800  continue
      print*,'Error 800, undefined LATS output convention'
      stop 800

 999  continue

      stop
      end
      
      subroutine read_data(a,name,ni,nj,id,it)
C         
C         routine to generate quasi realistic AMIP II fields
C         
      dimension a(ni,nj)
      character name*9

      pi=3.1459

C         
C         sfc field: psl = sea level pressure
C

      pscl=10.0
      pmin=950.0
      p0=1013.0

      i0=ni/2
      j0=nj/2
      rly=j0
      rlx=i0*0.5
C         
C         low eastward  in time
C
      i0=i0+(it-1)*3


      if(id.eq.0) then

        do i=1,ni
          do j=1,nj
            x=(i-i0)*1.0
            y=(j-j0)*1.0
            r=sqrt( x*x + y*y )
            p=p0-(p0-pmin)*exp(-r/pscl)
            a(i,j)=p*100.0
          end do
        end do

        return

      endif
C         
C         upper air field temperature
C 
C         850 mb
C         

      if(id.eq.1) then
        t0=250.0
        ay=10.0
        ax=30.0
        ishft=0
      endif

C         
C         500 mb
C 

      if(id.eq.2) then
        t0=235.0
        ay=20.0
        ax=20.0
        ishft=-5
      endif
C         
C         200 mb
C
      if(id.eq.3) then
        t0=210.0
        ay=30.0
        ax=10.0
        ishft=-10
      endif
C         
C         shift pattern eastward in time
C
      ishft=ishft-(it-1)*2

      do j=1,nj
        y=j0-(j-1)
        do i=1,ni
          angy=pi*0.5*(y/rly)
          angyx=angy*2
          x=i0-(i-1)+ishft
          angx=pi*0.5*(x/rlx)
          a(i,j)=t0 + ay*cos(angy) + ax*sin(angyx)*sin(angx)
        end do
      end do

      return
      end
