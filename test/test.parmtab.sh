#!/bin/sh

#f77 $1.f -L/usr/local/lats/LATS-1.1/lib -llats -L/usr/local/netcdf-2.3.2/libsrc -lnetcdf -o $1
f77 test.parmtab.f -L../ -llats -L/usr/local/lats/LATS-1.0 -lnetcdf
