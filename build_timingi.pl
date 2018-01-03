#!/usr/bin/perl -w
open(TFILE,">timingi.h") or die "Could not open timingi.h";
open(IFILE,"timingii.h") or die "Could not open timingii.h";
print TFILE "#include \"djb_timing.h\"\n";
print TFILE "#ifndef _PFPID_TIMING_\n";
print TFILE "#define _PFPID_TIMING_\n";
$i=3;
while (<IFILE>) {
    @defs = split;
    print TFILE "#define $defs[0] $i,$defs[1]\n";
    $i++;
}
print TFILE "#endif\n" ;
close TFILE;
close IFILE;
