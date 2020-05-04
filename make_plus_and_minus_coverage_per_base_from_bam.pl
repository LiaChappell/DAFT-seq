#!/usr/bin/perl -w
use strict;

unless (@ARGV >0) {
        &USAGE;
}


sub USAGE {

    die '

Usage: script.pl inbam genome_file_for_bedtools
 
Makes two coverage files (per base); wrapper for bedtools.
' . "\n";
}

my $inbam =$ARGV[0];
my $genome_file= $ARGV[1];

system "bedtools genomecov -d -strand + -ibam $inbam -g $genome_file > plus_strand_per_base_coverage.$inbam.txt";
system "bedtools genomecov -d -strand - -ibam $inbam -g $genome_file > minus_strand_per_base_coverage.$inbam.txt";

exit; 
