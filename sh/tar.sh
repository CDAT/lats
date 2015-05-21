#!/bin/sh

tfile="testlats.tar"
cd ..
rm -i $tfile
tar -cvf $tfile test/lats.sh
tar -uvf $tfile amip2.parm
tar -uvf $tfile test/lats4dum.txt
tar -uvf $tfile test/testlats.f
tar -uvf $tfile test/wgrib.c
tar -uvf $tfile test/latsout.grb.orig
tar -uvf $tfile sh/tar.sh
compress $tfile

exit
