#!/usr/bin/perl -w
open(TFILE,">fptimingi.h") or die "Could not open fptimingi.h";
open(IFILE,"fptimingii.h") or die "Could not open fptimingii.h";
print TFILE "#include \"djb_timing.h\"\n";
print TFILE "#ifndef _PFPID_TIMING_\n";
print TFILE "#define _PFPID_TIMING_\n";
$i=0;
while (<IFILE>) {
    @defs = split;
    print TFILE "#define $defs[0] $i,$defs[1]\n";
    $i++;
}
print TFILE "#endif\n" ;
close TFILE;
close IFILE;
