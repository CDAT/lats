      program testlats

      include "/d1/lats/lats.inc"

      parameter(ni=72,nj=46,nk=3,nv=2)

      character*20 center
      character*20 model
      character*9 var

      dimension var(nv),id_var(nv)
      double precision rlon(ni),rlat(nj),plev(nk)
      real ta(ni,nj,nk),psl(ni,nj)

      data center/'PCMDI'/
      data model/'lats'/
      data var/'psl','ta'/ 
      data plev/850.0,500.0,200.0/

      latsconv=LATS_PCMDI
      latsconv=LATS_GRIB_ONLY
      latsconv=LATS_GRADS_GRIB
      latsconv=LATS_COARDS

      latscal=LATS_NOLEAP
      latscal=LATS_CLIMLEAP
      latscal=LATS_CLIM
      latscal=LATS_STANDARD

      latsfreq=LATS_HOURLY
      latsfreq=LATS_FIXED
      latsfreq=LATS_MONTHLY

      latsdelta=0
      latsdelta=1

      id_var(1)=0
      id_var(2)=0

      if(latsconv.eq.LATS_GRADS_GRIB.or.
     $     latsconv.eq.LATS_GRIB_ONLY) then
        id_fil = latscreate('latsout.test1',
     $       latsconv,
     $       latscal,
     $       latsfreq,latsdelta,center,
     $       model,'LATS GRIB test')
        print*,'LATS grib file id = ',id_fil
      else if(latsconv.eq.LATS_COARDS.or.latsconv.eq.LATS_PCMDI) then
        id_fil = latscreate('latsout.test1',
     $       latsconv,
     $       latscal,
     $       latsfreq,latsdelta,center,
     $       model,'LATS netcdf test')
        print*,'LATS netcdf file id = ',id_fil

      endif

      nmo=2
      iyr0=1978
      imo0=1
      ida0=1
      ihr0=0

cc      if(latsfreq.ne.LATS_FIXED) then
cc        call latsbasetime(id_fil,iyr0,imo0,ida0,ihr0)
cc      endif


      do i=1,ni
        rlon(i)=0.0+(i-1)*360.0/ni
      end do

      do j=1,nj
        rlat(j)=-90.0+(j-1)*180.0/(nj-1)
      end do

      id_grd=latsgrid("u54", LATS_LINEAR, ni, rlon, nj, 
     $     rlat)
      print*,'id_grd = ',id_grd

      id_lev=latsvertdim('pressure', 'plev', nk, plev)
      print*,'id_lev = ',id_lev

      id_var(1)=latsvar(id_fil,var(1),LATS_FLOAT,LATS_AVERAGE,id_grd,
     $     0, 'sfc variable')
      print*,' id_var(1) = ',id_var(1)

      do imo=1,nmo

        if(latscal.eq.LATS_CLIM) then
          iyr=0
        else
          iyr=1979
        endif

        ida=1
        ihr=0

C         
C         create sfc field and write
C
        call read_data(psl,var(1),ni,nj,0,imo)
        ierr=latswrite(id_fil,id_var(1),0,iyr,imo,ida,ihr,psl)
        print*,'---------------------- ',imo,' ierr = ',ierr
        if(ierr.eq.0) stop 'latswrite error - sfc'


        if(id_var(2).eq.0) then


        id_var(2)=latsvar(id_fil,var(2),
     $       LATS_FLOAT,LATS_AVERAGE,id_grd,
     $       id_lev, 'ua variable')

c        id_var(2)=latsvar(id_fil,'tas',
c     $       LATS_FLOAT,LATS_AVERAGE,id_grd,
c     $       0, 'ua variable')

         print*,'id_var(2) = ',id_var(2)

      endif

        do k=1,nk
          call read_data(ta(1,1,k),var(2),ni,nj,k,imo)
          ierr=latswrite(id_fil,id_var(2),plev(k),iyr,imo,ida,ihr,
     $       ta(1,1,k))
          if(ierr.eq.0) stop 'latswrite error - ua'
        end do
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

C666666666666666666666666666666666666666666666666666666         
C         
C         STEP #6 - close the file
C         
C         1) for the GRIB, this creates the GrADS .ctl file for
C         futher processing by cdunif, GrADS or VCS
C
C666666666666666666666666666666666666666666666666666666         

      ierr = latsclose(id_fil)
      stop
      end
      
      subroutine read_data(a,name,ni,nj,id,it)
C         
C         routine to generate quasi realistic AMIP II fields
C         
      dimension a(ni,nj)


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
