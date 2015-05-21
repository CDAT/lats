#!/bin/sh

f77 testlats.f -L.. -llats -L/usr/local/lib -lnetcdf -o testlats

exit
