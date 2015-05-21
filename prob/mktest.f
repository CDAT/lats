C f77  -C -e  -g $1.f -L/usr/local/lats/lib -llats -L$PCMDI/lib -lezget -lcdms -ldrs -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -L/usr/local/HDF4.0r2/lib -lmfhdf -ldf -ljpeg -lz -o /pcmdi/doutriau/bin/$1
      PROGRAM testmkamip
      INTEGER nlon,nlat,ntime,len
      PARAMETER(nlon=360,nlat=180,ntime=12)
      REAL sst(nlon*nlat*ntime),wtfrac(nlon*nlat*ntime)
      CHARACTER*80 file,sftfileo
      CHARACTER*64 var,sftnameo
       INCLUDE '/usr/local/lats/include/lats.inc'
      CALL initget
      file='/pcmdi/staff/longterm/fiorino/sic/orig/amip2.bcs.sic.ctl'
      var='sica'
      sftfileo = '/pcmdi/staff/longterm/doutriau/ldseamsk/'//
     $     'amipII/owned/land_sea_mask_ukmo.nc'
      sftnameo='sftlf'
      CALL defvar (2,var,file)
        call defdim(2, 1, 'longitude', 'width', 'range', 
     &                        0.0, 0.0, 360.0)
        call defdim(2, 2, 'latitude', 'cosine', 'range',
     &                        0.0, 0.0, 0.0)
        call defdimi(2, 3, 'time', 'unit', 1, 1)
        call defdim(3, 1, 'longitude', 'width', 'as saved',
     &                        0.0, 0.0, 360.0)
        call defdim(3, 2, 'latitude', 'cosine', 'as saved',
     &                        0.0, 0.0, 0.0)
        call defdimi(3, 3, 'time', 'unit', 1, 1)

        call defvar(3, sftnameo, sftfileo)

        call defregrd(2, 'to', 3, 'area-weighted', 
     &                        0, 0.0, 0.0, 0, 0.0, 0.0)

        len = nlon*nlat
        call defmisc('data size', 'integer', len)

        print*, 'Entering EzGet subroutine getvdata'

        call getvdata(2, wtfrac, sst)

        print*, 'Successfully returned from EzGet subroutine getvdata'

      WRITE(9,*) (sst(i),i=1,100)
      CALL closeget

      i=latscreate('test.nc',
     $     LATS_GRADS_GRIB,LATS_STANDARD,LATS_FIXED,1,
     &     'center','model','com')
      END
