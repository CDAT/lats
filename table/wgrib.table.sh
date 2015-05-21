#!/bin/sh

l.pl $1.parms.txt wgrib | sort +0 -n -t\: > $1.wgrib.table

exit