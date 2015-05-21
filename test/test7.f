c f77 -e -g test7.f -L/pcmdi/ktaylor/pcmdi/util/ezget -lezgetdebug -L /usr/local/grads/lats -llats -L$PCMDI/lib -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o test7

c 3456789012345678901234567890123456789012345678901234567890123456789012


      parameter (nlon=6, nlat=3)
      implicit none
      integer i, j, nlon, nlat, igridtyp, icalend, ierr1, ifiletyp,
     &    igrid, istat, ifreq, inctime, idfile, ilevel, iyr, 
     &    imon, iday, ihour
      real omax(nlon,nlat), cmax(nlon,nlat) 
      double precision alev, alons(nlon), alats(nlat)
      character*80 varcomm, outfile
      character*16 varname, model
      character*120 filecomm
      character*16 gridtype, center 
      character*256 parmtabl
      integer idvara(2)
      integer latsgrid, latsparmtab, latscreate, latsvar, latswrite, 
     &     latsclose

      include '/usr/local/include/lats.inc'

c       ifiletyp = LATS_GRADS_GRIB, LATS_PCMDI, LATS_COARDS
      ifiletyp = LATS_GRADS_GRIB
      ifiletyp = LATS_GRIB
c      ifiletyp = LATS_COARDS

c       icalend = LATS_STANDARD, LATS_JULIAN, LATS_CLIM
      icalend = LATS_STANDARD

c       istat = LATS_INSTANT, LATS_AVERAGE
      istat = LATS_AVERAGE

      outfile =  'test7.out' 

    
      do 100 i=1,nlon
          alons(i) = (i-1)*60
  100 continue

      do 200 j=1,nlat
        alats(j) = (j-2)*60
        do 150 i=1,nlon
          omax(i,j) = 0.0
          cmax(i,j) = 0.0
  150   continue
  200 continue

      igridtyp = LATS_LINEAR
      gridtype = 'linear'
   
      print*, '    Defining grid with LATS ....'
  
      igrid = latsgrid(gridtype, igridtyp, nlon, alons, nlat, alats)

      if (igrid .eq. 0) then
        print*, 'Error in defining grid with latsgrid'
      else
        print*, 'Grid defined successfully by LATS'
      endif

      parmtabl = '/pcmdi/ktaylor/pcmdi/util/ketgrib.parms'

c      ierr1 = latsparmtab(parmtabl)

c      if (ierr1 .eq. 0) then
c        print*, 'Error -11 '
c      endif
                
      ifreq = LATS_FIXED
ccc      ifreq = LATS_MONTHLY
      inctime = 1
      model = 'test'
      center = 'pcmdi'
      filecomm = 'source 1'
         
      print*, ' '
      print*, '     Creating file with LATS:  ', outfile

      idfile = latscreate(outfile, ifiletyp, icalend, 
     &       ifreq, inctime, center, model, filecomm)

        if (idfile .eq. 0) then
          print*, 'Error in creating file in LATS '
        endif

        print*, '    Defining a variable with LATS: ', varname

      ilevel = 0
      varname = 'sstbcs'
      varcomm = 'title 1'

      idvara(1) = latsvar(idfile, varname, LATS_FLOAT, istat, igrid,
     &       ilevel,  varcomm)

        if (idvara(1) .eq. 0) then
          print*, 'Error in defining variable in LATS'
        endif


      varname = 'sst'
      varcomm = 'title 2'

c      idvara(2) = latsvar(idfile, varname, LATS_FLOAT, istat, igrid,
c     &       ilevel,  varcomm)

c        if (idvara(2) .eq. 0) then
c          print*, 'Error in defining variable in LATS'
c        endif
c
      alev = 0.0
      iyr = 1950
      imon = 1
      iday = 0
      ihour = 0

      ierr1 = latswrite(idfile, idvara(1), alev, iyr,
     &            imon, iday, ihour, omax)

        if (ierr1 .eq. 0) then
          print*, 'Error in writing field '
        endif

c      ierr1 = latswrite(idfile, idvara(2), alev, iyr,
c     &            imon, iday, ihour, cmax)

c        if (ierr1 .eq. 0) then
c          print*, 'Error in writing field '
c        endif

      ierr1 = latsclose(idfile)
        if (ierr1 .eq. 0) then
            print*, 'Error in closing file in LATS'
        endif

      call exit(1)

      end





