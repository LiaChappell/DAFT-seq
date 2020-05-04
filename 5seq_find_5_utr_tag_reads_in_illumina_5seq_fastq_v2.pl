#!/usr/bin/perl -w

use strict;

unless (@ARGV >0) {
        &USAGE;
}


sub USAGE {

    die '

Usage: script.pl blatout.psl 

Takes a blat-output, and retrieves the names of reads that are hits.

**Warning: psl-files should have no header**

' . "\n";
}


my $blat = shift;


# read in BLAT into a hash of hashes, with hit-name as key, and full line after
# merge overlapping hits to the longest

open (IN, "<$blat")|| die  "Where is the blat?\n";


# reading in the hits and filter bad ones

while (<IN>) {
    chomp;

    #print "$_\n";
    my @arr = split(/\s+/, $_);

    if ($_=~/^\d+/) { 
	if ($arr[0]> 11 and $arr[1] < 2 and $arr[5] < 2 and $arr[7] < 2 ){ 
		#print "$_\n";
		print "$arr[9]\n";
		}
	}
}


