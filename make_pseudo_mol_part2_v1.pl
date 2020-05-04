#!/usr/bin/perl -w
use strict;
use Getopt::Long;

## File inputs
my $outpref; # prefix of coverage files, e.g. 3D7_v3
my $thr;

GetOptions
(
	"p|pref:s"					=> \$outpref,
	"t|thr:s"					=> \$thr,
);

if (!defined $outpref ||!defined $thr ){
	print_usage();
	exit;
}


## Yell which inputs and threshold are in use
open (LOG, ">$outpref.thr$thr.log.txt") or die "$!";
print LOG "Inputs to script:\noutpref=ARGV[0]; # prefix of coverage files, e.g. 3D7_v3\nthr= ARGV[1]; # read coverage threshold, e.g. 5,10,15 .. 100\ngenes=ARGV[2]; #just the genes you want UTRs for, probably protein-coding genes\nall_feat=ARGV[3]; #include other annotated elements, such as ncRNAs etc\n\n";
print LOG "Starting analysis for a threshold of $thr reads\n";

## Get coverage at different ranges of depth- running this script once for each threshold
open (COV_P, "<$outpref.coverage.all.plus.strand.bed") or die "$!";
open (COV_P_GTET, ">$outpref.tmp.2.coverage.gtet.$thr.plus.strand.bed") or die "$!";
while (<COV_P>){
	my @covp= split /\s+/;
	if ($covp[3] > $thr){
		print COV_P_GTET "$_";
	}
}
close COV_P;
close COV_P_GTET;

open (COV_M, "<$outpref.coverage.all.minus.strand.bed") or die "$!"; 
open (COV_M_GTET, ">$outpref.tmp.2.coverage.gtet.$thr.minus.strand.bed") or die "$!";
while (<COV_M>){
	my @covm= split /\s+/;
	if ($covm[3] > $thr){
		print COV_M_GTET "$_";
	}
}
close COV_M;
close COV_M_GTET;

## Merge entries to give blocks of continuous coverage - "pseudo-molecules"
# Only positions which are next to each other- no gaps
# Merged bedfile has three columns: chr, feature start, feature end
system "bedtools merge -i $outpref.tmp.2.coverage.gtet.$thr.plus.strand.bed > $outpref.tmp.3.coverage.gtet.$thr.merged.plus.strand.bed";
system "bedtools merge -i $outpref.tmp.2.coverage.gtet.$thr.minus.strand.bed > $outpref.tmp.3.coverage.gtet.$thr.merged.minus.strand.bed";
print LOG "Merged coverage to make pseudo-molecules at read depth threshold of $thr\n";


## Turn temporary bed files into GFF files
open (BED_P, "<$outpref.tmp.3.coverage.gtet.$thr.merged.plus.strand.bed") or die "$!"; 
open (OUT_P, ">$outpref.coverage.gtet.$thr.merged.plus.strand.gff") or die "$!";
while(<BED_P>){
	chomp;
	my ($chr, $start_bed, $end)= split/\t/;
	my $start = ($start_bed +1); #bodge bed zero start to gff one start format
	print OUT_P "$chr\tPseudo_molecule_maker\tCoverage\t$start\t$end\t.\t+\t.\tParent=unknown\n";
}
close BED_P;
close OUT_P;

	
open (BED_M, "<$outpref.tmp.3.coverage.gtet.$thr.merged.minus.strand.bed") or die "$!";
open (OUT_M, ">$outpref.coverage.gtet.$thr.merged.minus.strand.gff") or die "$!";
while(<BED_M>){
	chomp;
	my ($chr, $start_bed, $end)= split/\t/;
	my $start = ($start_bed +1); #bodge bed zero start to gff one start format
	print OUT_M "$chr\tPseudo_molecule_maker\tCoverage\t$start\t$end\t.\t-\t.\tParent=unknown\n";
}
close BED_M;
close OUT_M;



##Tidy clutter
system "rm $outpref.tmp.*.$thr.*";
print LOG "Removed temporary files\n";

## Done with this threshold
print LOG "Done with munging for read threshold $thr\n";


sub print_usage
{
	print <<USAGE;
	
	Input files: 
        "p|pref:s"                   => Prefix of output from previous script- e.g. 3D7 

	Parameters: 
        "t|thr:s"                    => Threshold for coverage blocks **required**


USAGE
}
