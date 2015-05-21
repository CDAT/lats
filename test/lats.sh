#!/bin/sh

#
#	define the system variables
#
#
LATSVER="1.4"
BDIR="/home1/grads/src160"
BDIR="/usr/local"
NCLIB="$BDIR/lib"
NCLIB="$BDIR/lats/netcdf/lib"
NCLIB=""
LATSLIB="$BDIR/lats/lats-1.4/"
LATSLIB="$BDIR/lats"
LATSINC="$BDIR/lats/lats.inc"
F77="f77"
CC="cc"

#
#	typhoon.llnl.gov (sunos4)
#

LATSVER=""
BDIR="/d1"
NCLIB="/usr/local/lib"
LATSLIB="$BDIR/lats"
LATSLIB="$BDIR/lats"
LATSINC="$BDIR/lats/lats.inc"

WGRIB="$BDIR/test"
GRADS="/usr/local/grads/bin"
CUDUMP=""

ofile="latsout"
format="netcdf"
format="grib"

#
#	make the test routine by a sed of the template code testlats.f.mak
#	this looks weird but it works...
#

if [ $1 = "makecode" ] ; then
	
  echo "#!/bin/sh " > junk.sh
  echo "sed 's,XXXXX,$LATSINC,g' testlats.f.mak > testlats.f" >> junk.sh
  chmod a+x junk.sh
  junk.sh
  rm junk.sh

#
#	compile the routine
#
elif [ $1 = "c" ] ; then

  if [ $NCLIB ] ; then
    create="$F77 testlats.f -o testlats -L$LATSLIB -llats -L$NCLIB -lnetcdf"
  else
    create="$F77 testlats.f -o testlats -llats -L$LATSLIB"
  fi
  echo  "lllllllllllllllllllllllllllllll"
  echo " "
  echo "the command to create the LATS appliction testlats:"
  echo " "
  echo $create
  echo  " "
  echo  "lllllllllllllllllllllllllllllll"

  
  $create

#
#	run it
#
elif [ $1 = "makewgrib" ] ; then

  $CC wgrib.c -o wgrib -lm
#
#	run it
#
elif [ $1 = "run" ] ; then

  testlats

#
#	use wgrib to scan it
#
elif [ $1 = "wgrib" -a $format = "grib" ] ; then

  $WGRIB/wgrib "$ofile.grb"

#
#	create the gribmap for cdunif/grads/vcs
#
elif [ $1 = "gribmap" ] ; then

  $GRIBMAP/gribmap -v -i "$ofile.ctl"

#
#	run grads
#
elif [ $1 = "grads" ] ; then

  $GRADS/grads -lc "run $ofile.gs"

#
#	dump using cudump
#
elif [ $1 = "cudump" -a $format = "grib" ] ; then

  $CUDUMP/cudump $ofile.ctl

elif [ $1 = "cudump" -a $format = "netcdf" ] ; then

  $CUDUMP/cudump $ofile.nc

fi


exit
