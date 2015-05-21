#!/usr/local/bin/perl

$perldir=$ENV{"WXMAP_PERL_DIR"};
require("$perldir/lib/mf.pl");
$wxdir=$ENV{"WXMAP_PRC_WXMAP_DIR"};
require("$wxdir/wxmap.env.pl");
$webmast=wxmap_master();

$file=$ARGV[0];
$opt=$ARGV[1]; 

$narg=$#ARGV+1;

if($narg >= 1) {
  $file=$ARGV[0];
  $opt=$ARGV[1]; 
} else {

  print "\n";
  print "$0 processing :\n\n";
  print "     file : create_data | input file\n";
  print "      opt : wgrib | amip2 | obs.var.htm\n";
  print "  Try again\n\n";
  exit;

}

if($file eq "create_date") {
  print_create_date();
  exit;
}

#print "$file\n";
open(FILE,$file) || die "Can't open $file\n" ;

while(<FILE>) {
    $card=$_;
    @c=split('\|',$card);
    @t=split(',',$c[8]);
    if ($t[0] != '' || $opt eq "obs.var.htm") { 
      $card="$card";
      push(@n,$card);
    }
#    print "CARD:  $c[0] $c[1] $c[9] \n";
}

@n=sort by_varnum @n if($opt eq "sort_by_varnum") ;
@n=sort by_varname @n if($opt =~ "amip2" || $opt eq "obs.var.htm") ;

var_table() if($opt eq "obs.var.htm");

if($opt eq "amip2") {

foreach $nn (@n) {
  print "$nn";
}

}

exit;

if($opt eq "wgrib") { print "-1:100:-1:-1\n"; }
@n=sort(@n);
for( $i = 0; $i <= $#n; $i++ ) {

   if($opt eq "wgrib") {
     @c=split('\|',$n[$i]);
     printf("%s:%s:%s [%s]\n",$c[1],$c[0],$c[2],$c[3]);
   } else {
     printf("%s",$n[$i]);
   }

  }

exit;


sub var_table {

$hfile="obs.mo.var.htm";
$htitle="Variables in /pcmdi/obs/mo";

$bdir="/pcmdi/obs/mo";

$curdir=`pwd`;
chomp($curdir);

chdir($bdir);

$i=0;
while ( <*> ) {
  
    $vd[$i] =$_ ;
    print "$vd[$i]\n" if(-d $vd[$i]);
    $i++;
}


chdir($curdir);


#
#  create the HTML
#

open(H,"> $hfile") || die "Can't open $hfile\n" ;

$width{'name'}=75;
$width{'desc'}=325;
$width{'units'}=100;

$hbody="
<html>
<head>
<title>$htitle</title>
</head>
<h1><i><font face=\"arial\">$htitle</font></i></h1>
<br>
<br>
<table border=1 cellpadding=5 cellspacing=2>
<tr>
<th width=$width{'name'} align=center valign=center>
Name
</th>
<th width=$width{'desc'} align=left valign=center>
Description
</th>
<th width=$width{'units'} align=left valign=center>
Units
</th>
</tr>
";

print H $hbody;

foreach $nn (@n) {
  @c=split('\|',$nn);

  foreach $vvd (@vd) {
    $vvd;
    @cc=split(' ',$c[0]);

    if($cc[0] eq $vvd) {

      printf("%s:%s:%s [%s]\n",$c[1],$c[0],$c[2],$c[3]);


  $hrow="
<tr>

<td width=$width{'name'} align=center valign=center>
<b>$c[0]</b>
</td>

<td width=$width{'desc'} align=left valign=center>
$c[2]
</td>

<td width=$width{'units'} align=left valign=center>
$c[3]
</td>

</tr>
";

  print H $hrow;

  last;
  
    }

  }

}
$htail="
</tr>
</table>
<br>
<br>
$webmast
</body>
</html>
";

print H $htail;
  
close(H);

}

sub print_create_date {

$curtime=dtg("timesec");
$project="AMIP II"; 

$tt="
#---------------------------------------------------------------------------------------------------
#
#  This LATS parmeter file for $project was created:  $curtime
#
#  by Mike Fiorino, PCMDI fiorino\@pcmdi.llnl.gov
#
#---------------------------------------------------------------------------------------------------

";

print $tt;
}

#
# compare routine for sorting by parm number element #1 in array
#

sub by_varnum {
  $ii=1;
  my($aa,$bb);
  @aa=split(/\|/,$a);
  @bb=split(/\|/,$b);
  ($aa[$ii] <=> $bb[$ii]) || ($aa[0] cmp $bb[0]) ;

}

#
# compare routine for sorting by parm number element #1 in array
#

sub by_varname {
  my($aa,$bb);
  @aa=split(/\|/,$a);
  @bb=split(/\|/,$b);
  ($aa[0] cmp $bb[0]) ;

}
