#!/bin/sh

ver="0.4"
h=`echo $HOST | awk -F. '{print $1}'`

DEBUG=""

if [ $# -eq 1 ] ; then
  DEBUG="-g"
fi

if [ $h = "tenki" ] ; then
  cc $DEBUG latsls.c latsint.o -o latsls -DLATSLS_VERSION=\"${ver}\" -lm
else
  acc $DEBUG latsls.c latsint.o -o latsls -DLATSLS_VERSION=\"${ver}\" -lm
fi
