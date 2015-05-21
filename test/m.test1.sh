#!/bin/sh

f77 -g test1.f -o test1 -L.. -llats -L/usr/local/lib -lnetcdf

exit
