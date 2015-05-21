      program testsvert

CCC       f77  $1.f -L./ -llats -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o $1

      include '/usr/local/lats/LATS-1.1/include/lats.inc'
      
      integer latsvertdim, latscreate, latsvar, latswrite,latsclose
      integer latsgrid
      integer grid_72x96
      
      integer latsfile, ierr
      integer dpuv, dpt, dpall
      integer ibyr, ibmon, ibdat, ibhr

      integer ua, ta
      
      integer mgrid
      parameter (mgrid=72*96)

      double precision lons(96) /
     &     0.,      3.75,    7.5,    11.25,   15.,     18.75,
     &     22.5,    26.25,   30.,
     $     33.75,   37.5,    41.25,   45.,     48.75,   52.5,
     $     56.25,   60.,     63.75,
     &     67.5,    71.25,   75.,     78.75,   82.5,    86.25,
     $     90.,     93.75,   97.5, 
     &     101.25,  105.,    108.75,  112.5,   116.25,  120.,
     $     123.75,  127.5,   131.25,
     &     135.,    138.75,  142.5,   146.25,  150.,    153.75,  
     $     157.5,   161.25,  165.,  
     &     168.75,  172.5,   176.25,  180.,    183.75,  187.5,
     $     191.25,  195.,    198.75,
     &     202.5,   206.25,  210.,    213.75,  217.5,   221.25,
     $     225.,    228.75,  232.5, 
     &     236.25,  240.,    243.75,  247.5,   251.25,  255.,    
     $     258.75,  262.5,   266.25,
     &     270.,    273.75,  277.5,   281.25,  285.,    288.75,
     $     292.5,   296.25,  300.,  
     &     303.75,  307.5,   311.25,  315.,    318.75,  322.5,
     $     326.25,  330.,    333.75,
     &     337.5,   341.25,  345.,    348.75,  352.5,   356.25 /
      
      double precision lats(72) /
     &     88.75,  86.25,  83.75,  81.25,  78.75,  76.25,  73.75,  
     $     71.25, 68.75,  66.25,  63.75,
     &     61.25,  58.75,  56.25,  53.75,  51.25,  48.75,  46.25,
     $     43.75,  41.25,  38.75,
     &     36.25,  33.75,  31.25,  28.75,  26.25,  23.75,  21.25,
     $     18.75,  16.25,  13.75,
     &     11.25,   8.75,   6.25,   3.75,   1.25,  -1.25,  -3.75,
     $     -6.25,  -8.75, -11.25,
     &     -13.75, -16.25, -18.75, -21.25, -23.75, -26.25, -28.75,
     $     -31.25, -33.75, -36.25,
     &     -38.75, -41.25, -43.75, -46.25, -48.75, -51.25, -53.75,
     $     -56.25, -58.75, -61.25,
     &     -63.75, -66.25, -68.75, -71.25, -73.75, -76.25, -78.75,
     $     -81.25, -83.75, -86.25,
     &     -88.75 /

c     Define some fake data

      real x850(mgrid) /mgrid*850.0/
      real x500(mgrid) /mgrid*500.0/
      real x200(mgrid) /mgrid*200.0/
      real x50(mgrid) /mgrid*50.0/

      double precision dpuvvals(3) / 850.0, 200.0, 50.0 /
      double precision dpallvals(4) / 850.0, 500.0, 200.0, 50.0 /
      double precision dptvals(3)  / 850.0, 500.0, 50.0 /

      logical doerror

      doerror=.false.
      doerror=.true.

c     Define the grid(s)
      
      grid_72x96 = latsgrid('72x96', LATS_LINEAR, 96, lons, 72, lats)
      if (grid_72x96.eq.0) stop 'Grid def error'

c     Define the vertical levels

      if(doerror) then
        dpuv = latsvertdim('dpuv', 'plev', 3, dpuvvals)
        dpt  = latsvertdim('dpt', 'plev', 3, dptvals)
      else
        dpall = latsvertdim('dpall', 'plev', 4, dpallvals)
        dpuv=dpall
        dpt=dpall
      endif

      latsfile = latscreate('twolevs',
     $     LATS_GRADS_GRIB,LATS_STANDARD,LATS_HOURLY,6,'ukmo',
     &  'hadAM3, 19 levels','Test')
      if (latsfile.eq.0) stop 'Error 6-hourly (B stream) file'

      ua = latsvar(latsfile, 'ua',
     $     LATS_FLOAT, LATS_INSTANT, 
     $     grid_72x96, dpuv, 'u compnt of wind on pressure levels   ')

C         
C
C *** THIS SHOULD GENERATE AN ERROR ! ***
C         mf 970822 - done         
C

      ta = latsvar(latsfile, 'ta',
     $     LATS_FLOAT, LATS_INSTANT, 
     $     grid_72x96, dpt,
     $     't on pressure levs u grid. use macr   ')

      ibyr = 1978
      ibmon = 12
      ibdat = 1
      ibhr = 0

c     Write ua at 850, 200, 50 hPa
      
      ierr = latswrite(latsfile, ua,
     $     850.0D0, ibyr, ibmon, ibdat, ibhr, x850)
      ierr = latswrite(latsfile, ua,
     $     200.0D0, ibyr, ibmon, ibdat, ibhr, x200)
      ierr = latswrite(latsfile, ua,
     $     50.0D0, ibyr, ibmon, ibdat, ibhr, x50)
      
c     Write ta at 850, 500, 50hPa

      ierr = latswrite(latsfile, ta, 
     $     850.0D0, ibyr, ibmon, ibdat, ibhr, x850)
      ierr = latswrite(latsfile, ta, 
     $     500.0D0, ibyr, ibmon, ibdat, ibhr, x500)
      ierr = latswrite(latsfile, ta, 
     $     50.0D0, ibyr, ibmon, ibdat, ibhr, x50)

      ierr = latsclose(latsfile)

      end

