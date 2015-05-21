#!/usr/local/bin/perl
# /* -*-Mode: perl;-*-
#
# Generate latsparm.h from a parameter file
#
# Usage: genlatsparm parameter_file
#
# $Id: genlatsparm.pl,v 1.5 1996/10/22 19:04:59 fiorino Exp $
#
# $Log: genlatsparm.pl,v $
# Revision 1.5  1996/10/22  19:04:59  fiorino
# latsgrib bug in .ctl creator
#
# Revision 1.4  1996/10/10 23:15:41  drach
# - lats_create filetype changed to convention, with options LATS_PCMDI,
#   LATS_GRADS_GRIB, and LATS_COARDS.
# - monthly data defaults to 16-bit compression
# - LATS_MONTHLY_TABLE_COMP option added to override 16-bit compression
# - AMIP II standard parameter file
# - parameter file incorporates GRIB center and subcenter
# - if time delta is positive, check that (new_time - old_time)=integer*delta
#
# Revision 1.3  1996/08/20 18:34:04  drach
# - lats_create has a new argument: calendar
# - lats_grid: longitude, latitude dimension vectors are now double
#   precision (double, C).
# - lats_vert_dim: vector of levels is now double precision (double,
#   C). lats_vert_dim need not be called for single-value/surface
#   dimensions, if defined in the parameter table. Multi-valued vertical
#   dimensions, such as pressure levels, must be defined with
#   lats_vert_dim.
# - lats_var: set level ID to 0 for implicitly defined surface
#   dimension.
# - lats_write: level value is double precision (double, C).
# - lats_parmtab: replaces routine lats_vartab.
# - FORTRAN latserropt added: allows program to set error handling
#   options.
# - The parameter file format changed.
#
# Revision 1.2  1996/06/27 01:01:08  drach
# - Removed time stats
# - changed min/max to range
# - handled empty tables properly
#
# Revision 1.1  1996/06/12 00:39:07  drach
# - Initial version
#

if($#ARGV != 0){
    &usage;
    exit 0;
}

open(IN, $ARGV[0]) || die "Can't open parameter file $ARGV[0]: $!\n";
if (-s "latsparm.h"){
    print "Save current latsparm.h? [y/n]:";
    $answer = scalar(<STDIN>);
    $save = ($answer !~ /^n?$/i);
    if($save){
	rename("latsparm.h","latsparm.h.$$") || die "Can't rename latsparm.h\n";
	print "Saved latsparm.h as latsparm.h.$$\n";
    }
}
open(OUT,">latsparm.h") || die "Can't create latsparm.h: $!\n";

				# Print the initial output lines
print OUT '/* -*-Mode: C;-*-',"\n";
print OUT ' * Module:      LATS default parameter table',"\n";
print OUT ' *',"\n";
print OUT ' * Copyright:	1996, Regents of the University of California',"\n";
print OUT ' *		This software may not be distributed to others without',"\n";
print OUT ' *		permission of the author.',"\n";
print OUT ' *',"\n";
print OUT ' * Author:      Bob Drach, Lawrence Livermore National Laboratory',"\n";
print OUT ' *              drach@llnl.gov',"\n";
print OUT ' *',"\n";

($sec,$min,$hour,$mday,$mon,$year) = (gmtime);
$mon++;
$login = getlogin;

printf OUT " * Generated:   $year-$mon-$mday $hour:$min:$sec ($login) from $ARGV[0] by genlatsparm.pl\n";
print OUT ' *',"\n";
print OUT ' * Version:     $Id: genlatsparm.pl,v 1.5 1996/10/22 19:04:59 fiorino Exp $',"\n";
print OUT ' */',"\n";
print OUT '#ifndef _LATSPARM_H',"\n";
print OUT '#define _LATSPARM_H',"\n";
print OUT '#include "latsint.h"',"\n";


				# Process each line of the input parameter file
$line = 0;
$mode = "start";
$firstblock = 1;
LINE:
while (<IN>){
    $line++;
    chop;
				# Look for comments and/or mode switches
    if (/^#/) {
	if (/^#!\s*(variable|vert|center|qc)/i){
	    &endblock($struct,$arrname,$macro) unless $firstblock;
	    $mode = $1;
	    $mode =~ tr/A-Z/a-z/;
	    if($mode eq "variable"){
		$struct = "latsParm";
		$arrname = "latsDefaultParms";
		$macro = "LATS_DEFAULT_NPARMS";
		$default_line = "{\"_default_\"}";
		&startblock($struct,$arrname);
		next LINE;
	    }
	    elsif($mode eq "vert"){
		$struct = "latsVertType";
		$arrname = "latsDefaultVerts";
		$macro = "LATS_DEFAULT_NVERTS";
		$default_line = "{\"_default_\"}";
		&startblock($struct,$arrname);
		next LINE;
	    }
	    elsif($mode eq "center"){
		$struct = "latsCenter";
		$arrname = "latsDefaultCenters";
		$macro = "LATS_DEFAULT_NCENTERS";
		$default_line = "{\"_default_\"}";
		&startblock($struct,$arrname);
		next LINE;
	    }
	    elsif($mode eq "qc"){
		$struct = "latsParmQC";
		$arrname = "latsDefaultQCs";
		$macro = "LATS_DEFAULT_NQCS";
		$default_line = "{\"_default_\"}";
		&startblock($struct,$arrname);
		next LINE;
	    }
	}
    }
				# Ignore blank line
    elsif (/^\s*$/) {
    }
				# Variables
    elsif ($mode eq "variable"){
	($name, $id, $title, $units, $datatype, $leveltype, $decimal_scale_factor, $precision, $comments_1, $comments_2) = split(/\s*\|\s*/);
	$name =~ tr/a-zA-Z0-9_//cd; # Remove leading whitespace missed by split

				# Map to enums
	$datatype =~ tr/A-Z/a-z/;
	if($datatype eq "float"){
	    $dtype = "LATS_FLOAT";
	}
	elsif ($datatype eq "int"){
	    $dtype = "LATS_INT";
	}
	else {
	    &badvalue("datatype",$datatype,$line);
	}

	if($leveltype eq ""){
	    $levelset = 0;
	}
	else{
	    $levelset = 1;
	}

				# Output the struct entry
	&startrec;
	&outstr($name);
	&outstr($title);
	&outstr($units);
	&outsc($id);
	&outsc($decimal_scale_factor);
	&outsc($precision);
	&outstr($leveltype);
	&outsc($levelset);
	&outsc(0);		# init pointer to verttype
	&outsc($dtype);
	&endrec;

				# Save the level type to compare with the vertical types table
	if($leveltype ne ""){
	    $surfs{$name} = $leveltype;
	}
    }
				# Vertical dimension types
    elsif ($mode eq "vert"){
	($level , $description,$units,$verticality,$positive,$grib_id,$grib_p1,$grib_p2,$grib_p3) = split(/\s*\|\s*/);
	$level =~ tr/a-zA-Z0-9_//cd; # Remove leading whitespace missed by split
				# Map to enums
	$verticality =~ tr/A-Z/a-z/;
	$positive =~ tr/A-Z/a-z/;
	if($verticality eq "single"){
	    $vert = "LATS_SINGLE_LEVEL";
	}
	elsif($verticality eq "multi"){
	    $vert = "LATS_MULTI_LEVEL";
	}
	else {
	    &badvalue("verticality",$verticality,$line);
	}
	if($positive eq "up"){
	    $pos = "LATS_UP";
	}
	elsif($positive eq "down"){
	    $pos = "LATS_DOWN";
	}
	else{
	    &badvalue("positive",$positive,$line);
	}
				# Output the struct entry
	&startrec;
	&outstr($level);
	&outstr($description);
	&outstr($units);
	&outsc($pos);
	&outsc($vert);
	&outsc($grib_id);
	&outsc($grib_p1);
	&outsc($grib_p2);
	&outsc($grib_p3);
	&endrec;

				# Save the verticalities and levels to compare with the variable table
	$leveltypes{$level} = $vert;
    }
				# Centers
    elsif ($mode eq "center"){
	($center,$grib_id,$grib_center,$grib_subcenter) = split(/\s*\|\s*/);
	$center =~ tr/a-zA-Z0-9_\///cd; # Remove leading whitespace missed by split
				# Output the struct entry
	&startrec;
	&outstr($center);
	&outsc($grib_id);
	&outsc($grib_center);
	&outsc($grib_subcenter);
	&endrec;
    }
				# QC marks
    elsif ($mode eq "qc"){
	($variable,$level_type,$level,$mean,$std,$tolerance,$range,$rangetol) = split(/\s*\|\s*/);
	$variable =~ tr/a-zA-Z0-9_\///cd; # Remove leading whitespace missed by split
	if($level eq ""){
	    $level=0.0;
	}
				# Output the struct entry
	&startrec;
	&outstr($variable);
	&outstr($level_type);
	&outsc($level);
	&outsc($mean);
	&outsc($std);
	&outsc($tolerance);
	&outsc($range);
	&outsc($rangetol);
	&outsc(0.0);		# Calculated mean
	&outsc(0.0);		# Calculated max
	&outsc(0.0);		# Calculated min
	&outsc(0.0);		# Calculated std
	&endrec;
    }
}

				# Last output lines
&endblock($struct,$arrname,$macro) unless $firstblock;
print OUT "#endif\n";

close IN;
close OUT;

				# Check that all the surfaces are defined in the vertical types table
foreach $varname (sort keys(%surfs)){
    if(!defined($leveltypes{$surfs{$varname}})){
	print "Error: variable \'$varname\' has surface \'$surfs{$varname}\'\n";
	print "  which does not appear in the vertical types table: FIX THIS!!!\n";
    }
				# and check that none of the variable surfaces are multi-level
    elsif($leveltypes{$surfs{$varname}} eq "LATS_MULTI_LEVEL"){
	print "Error: variable \'$varname\' has surface \'$surfs{$varname}\'\n";
	print "  which is multi-level: leave the surface blank or set to a single-level vertical dimension!!!\n";
    }
}

exit 0;

				# Start printing a static initialization block
sub startblock {
    local($struct,$arrname) = @_;
    print OUT "static $struct $arrname\[\] = {\n";
    $firstblock = 0;
    $nrec = 0;
}
				# End the current initialization block
sub endblock {
    local($struct,$arrname,$macro) = @_;
    print OUT $default_line if $nrec == 0;
    print OUT "\n};\n#define $macro (sizeof $arrname / sizeof($struct))\n"
}


				# Start the record output
sub startrec {
    print OUT ",\n" unless $nrec == 0;
    print OUT "{";
    $nfield = 0;
}
				# End the record
sub endrec {
    print OUT "}";
    $nrec++;
}
				# Print a string (enclose arg in quotes)
sub outstr {
    local($string) = @_;
    print OUT "," if $nfield > 0;
    print OUT "\"$string\"";
    $nfield++;
}
				# Print a scalar (no quotes)
sub outsc {
    local($scalar) = @_;
    print OUT "," if $nfield > 0;
    print OUT $scalar;
    $nfield++;
}
				# Report a bad value of attribute 'name'
				# with 'value' at 'line' of input file
sub badvalue {
    local($name,$val,$line) = @_;
    
    print "Error at line $line: invalid $name \'$val\'\n";
    exit(1);
}

sub usage {
    print "Usage: genlatsparm.pl <ASCII parameter file>\n";
}
