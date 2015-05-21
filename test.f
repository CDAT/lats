	program test
c  Test LATS.  This program is derived from .../NRL/readwrite.F.
c
c		Curt Covey	January 1997
c
c  To compile on SGI:
c
c	f77 -col120 test.F -I/usr/local/include \
c                          -L/usr/local/grads/lats -llats \
c		           -L/usr/local/lib -lnetcdf 
        include "lats.inc"
	parameter(nyr=36,mm=nyr*12)
        parameter(im=144,jm=72,lm=18)
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

c	 72  points in y:  

	data lats /
     .  -88.10,  -85.64,  -83.16,  -80.68,  -78.20,  -75.72,  -73.24,  -70.75,  -68.27,   
     .  -65.79,  -63.31,  -60.83,  -58.34,  -55.86,  -53.38,  -50.89,  -48.41,  -45.93,  -43.45,   
     .  -40.96,  -38.48,  -36.00,  -33.52,  -31.03,  -28.55,  -26.07,  -23.59,  -21.10,  -18.62,   
     .  -16.14,  -13.65,  -11.17,  -8.689,  -6.207,  -3.724,  -1.241,   1.241,   3.724,   6.207,   
     .   8.689,   11.17,   13.65,   16.14,   18.62,   21.10,   23.59,   26.07,   28.55,   31.03,   
     .   33.52,   36.00,   38.48,   40.96,   43.45,   45.93,   48.41,   50.89,   53.38,   55.86,   
     .   58.34,   60.83,   63.31,   65.79,   68.27,   70.75,   73.24,   75.72,   78.20,   80.68,   
     .   83.16,   85.64,   88.10/ 

c	 144 points in x:
c	      60.0E, 62.5E, 65.0E, ..., 57.5E

	print *, "Setting longitude values . . ."
	lons(1) = 60.0
	do i = 2, im
	   lons(i) = lons(i-1) + 2.5
	enddo
	print *, "Final value is ", lons(im)

	print *, "Creating LATS grid, file and variable . . ."
	idgrid = latsgrid("NRL_grid", LATS_GAUSSIAN, im, lons, jm, lats)
        idfile = latscreate("nrl", LATS_GRIB, LATS_STANDARD, LATS_MONTHLY, 1,
     .		            "nrl", "Coupled Ocean-Atmosphere GCM", "prescribed seaice")
        idvarb = latsvar(idfile, "tas", LATS_FLOAT, LATS_AVERAGE, idgrid, 0,
     .                   "Precise definitions of near-surface will vary.")
	if (idvarb .ne. 0) print *, "Creations appear to be successful."

	print *, "Putting 0's into test array . . ."
        nmo=1
	do k = 1, nmo
	   do j = 1, jm
	      do i = 1, im
	         ts(i, j, k) = 0.
	      enddo
	   enddo
	enddo
	print *, "Sample value: ts(72, 36, 216) = ", ts(72, 36, 216)

	print *, "Writing to LATS file . . ."
	icount = 0
	do iyr = 1900, 1900
	   print *, " . . . for Year ", iyr, " . . ."
	   do imo = 1, nmo
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
	   print *, "Writing appears to have been successful."
	else
	   print *, "WARNING: Writing appears to have been unsuccessful."
	endif

9901    format(10e13.6)

	irc=latsclose(idfile)
        print*,'irc = ',irc

c*cc    stop
	end
