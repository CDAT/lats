#!/bin/sh
#
#	clean up the variable listing
#
file=$1"
#cat $1 | awk -F\| '{print $1}'

#'{printf("%8s |%3s| %-50s |%10s |%5s|%6s|%6s|%6s| \n ",$1,$2,substr($3,1,50),$4,$5,$6,$7,$8)}'

sort -t\| +1 $file

# | awk -F\| '{printf("%8s|%3s|%-50s|%10s | %5s | %6s | %6s | %6s|\n ",$1,$2,substr($3,1,50),$4,$5,$6,$7,$8)}'

