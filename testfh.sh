#!/bin/sh

f77 -I. testfh.f liblats.a -L./lib -lnetcdf

exit
