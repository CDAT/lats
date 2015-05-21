#!/bin/sh

#
#	copy files to the local directory for users to access
#

tdir="/usr/local/grads/lats/test"
lfiles="lats.sh latsout.ctl latsout.grb latsout.nc latsout.gs"
tfiles="testlats.f testlats.f.mak lats4dum.txt"
cp $lfiles $tdir
cp $tfiles $tdir

exit
