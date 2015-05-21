#!/usr/local/bin/perl

$perldir=$ENV{"WXMAP_PERL_DIR"};
require("$perldir/lib/mf.pl");

$file=$ARGV[0];
$opt=$ARGV[1]; 

$narg=$#ARGV+1;

if($narg >= 1) {
  $file=$ARGV[0];
  $opt=$ARGV[1] if($narg >= 2); 
  $opt2=$ARGV[2] if($narg >= 3); 
} else {

  print "\n";
  print "$0 processing :\n\n";
  print "      file : create_data | input file\n";
  print "     [opt] : sort_by_varname | sort_by_varnum | wgrib | amip2\n";
  print "    [opt2] : short\n";
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
    if ($t[0] != '') { 
      $card="$card";
      push(@n,$card);
    }
#    print "CARD:  $c[0] $c[1] $c[9] \n";
}

@n=sort by_varnum @n if($opt eq "sort_by_varnum") ;
@n=sort by_varname @n if($opt eq "amip2" || $opt eq "sort_by_varname" ) ;

foreach $nn (@n) {
  if($opt2 eq "short") {
    @nnn=split('\|',$nn);
    print "$nnn[0] $nnn[1]\n";
  } else {
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
