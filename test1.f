	program readwrite

c  Read surface air temperature time series from NRL and write to NetCDF.
c
c		Curt Covey	January 1997
c
c  To compile on SGI:
c
c	f77 -col120 readwrite.F -I/usr/local/grads/lats \
c                               -L/usr/local/grads/lats -llats \
c			        -L/usr/local/lib -lnetcdf 

	integer testmode	! yes/no (1/0) switch for running a quick test
	data    testmode /1/

        include "lats.inc"

	parameter(nyr=36,mm=nyr*12)
ccc        parameter(im=144,jm=72,lm=18)
        parameter(im=72,jm=45,lm=18)
        dimension pr(im,jm),ps(im,jm),taux(im,jm),tauy(im,jm),
     .            hf(im,jm),wf(im,jm),seaice(im,jm),
     .            snow(im,jm)
        dimension tzm(jm,lm),uzm(jm,lm),vzm(jm,lm),qzm(jm,lm)
 	dimension ts(im,jm,mm)

	dimension tas(im,jm)	! This variable will be written to LATS.

c*cc    dimension sd1(jm,lm),sd2(jm,lm)

c=================================================================
c
c	NRL Coupled GCM Experiments (36-yr run)
c
c=================================================================
c
c	pr	precipitation (cm/day)
c	ps	sea-level pressure (mb)
c	taux	surface zonal wind stress (N/m2)
c       tauy    surface meridional wind stress (N/m2)
c	hf	net surface heat flux (W/m2)
c	wf	net surface freshwater flux (pr-evaporation) (cm/day)
c	seaice	sea ice cover (1: sea ice; 0: no sea ice)
c		(Note: seaice temperatures are calculated interactively in the model,
c		 but seaice cover is specified from climatology)
c	snow	snow depth (cm)
c	tzm	zonal mean air temperature (degree K)
c	uzm	zonal mean zonal wind (m/s)
c	vzm	zonal mean meridional wind (m/s)
c	qzm	zonal mean specific humidity (g/kg)
c	ts	36-yr (432-mth) lowest-layer air temperature (degree K)
c	sd1	standard deviation of annual-mean zonal-mean temperature (degree K)
c	sd2	standard deviation of decadal-mean zonal-mean temperature (degree K)
c
c----------------------------------------------------------------------------------------------
 
c	Atmospheric model grid:

	double precision lats(jm), lons(im)

	lons(1) = 0.0
	do i = 2, im
          lons(i) = lons(i-1) + 5.0
        enddo
        lats(1)=-88.0
        do j=2,jm
          lats(j)=lats(j-1)+4.0
        end do
        

	print *, "Final value is ", lons(im),lats(jm)

	print *, "Creating LATS grid, file and variable . . ."
	idgrid = latsgrid("NRL_grid", LATS_GENERIC, im, lons, jm, lats)
ccc	idgrid = latsgrid("NRL_grid", LATS_GAUSSIAN, im, lons, jm, lats)
        idfile = latscreate("nrl", LATS_COARDS, LATS_STANDARD, LATS_MONTHLY, 1,
     .		            "nrl", "Coupled Ocean-Atmosphere GCM", "prescribed seaice")
        idvarb = latsvar(idfile, "tas", LATS_FLOAT, LATS_AVERAGE, idgrid, 0,
     .                   "Precise definitions of near-surface will vary.")
	if (idvarb .ne. 0) print *, "Creations appear to be successful."

	if (testmode .eq. 0) then
	   print *, "Reading atmospheric data . . ."
           open(9992,file='/pcmdi/aux2/triton1/TwoDFields/NRL/cmip_atm.dta')

c for DJF
	   read(9992,9901) pr,ps,taux,tauy,hf,wf,
     .		seaice,snow,tzm,uzm,vzm,qzm

c for JJA
           read(9992,9901) pr,ps,taux,tauy,hf,wf,
     .          seaice,snow,tzm,uzm,vzm,qzm

           read(9992,9901) ts
c*cc       read(9992,9901) sd1,sd2
	else
           print *, "Putting 0.9's into test array . . ."
ccc           do k = 1, mm
           do k = 1, 1
              do j = 1, jm
                 do i = 1, im
                    ts(i, j, k) = 0.9
                 enddo
              enddo
           enddo
	endif
	print *, "Sample value: ts(72, 36, 216) = ", ts(72, 36, 216)

	print *, "Writing to LATS file . . ."
	icount = 0
	if (testmode .eq. 0) then
	   iyrmax = nyr 
	else
	   iyrmax = 1
	endif
	do iyr = 1, iyrmax
	   print *, " . . . for Year ", iyr, " . . ."
ccc	   do imo = 1, 12
	   do imo = 1, 1
	      icount = icount + 1
	      do j = 1, jm
	         do i = 1, im
	            tas(im, jm) = ts(im, jm, icount)
		 enddo
	      enddo
	      idreturn = latswrite(idfile, idvarb, 0.0D00, iyr, imo, 1, 0, tas)
           enddo
	enddo
	print *, "Total count of months = ", icount
	if (idreturn .eq. 1) then
	   print *, "Writing appears to have been successful, but file isn't closed yet."
	else
	   print *, "WARNING: Writing appears to have been unsuccessful."
	endif

	idreturn = latsclose(idfile)
        if (idreturn .eq. 1) then
           print *, "File appears to have been closed successfully."
        else
           print *, "WARNING: Closing appears to have been unsuccessful."
           print *, "         Return code from latsclose = ", idreturn
        endif

9901    format(10e13.6)

c*cc    stop
	end
