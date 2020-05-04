#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $put_ipt_gff; #made by another script
my $us_os_win; #upstream opposite strand windows
my $five_utr_gff; #coverage-based 5'UTRs
my $genome_file; #bedtools


GetOptions
(
	"i|put_ipt:s" 		=>\$put_ipt_gff,
	"u|us_on_win:s"		=>\$us_os_win,
	"f|five_utr_gff:s" 	=>\$five_utr_gff,
	"g|genome_file:s" 	=>\$genome_file,

);

sub print_usage
{
	print <<USAGE;
		
	"i|put_ipt:s" 		=>\$put_ipt_gff,
        "u|us_on_win:s"         =>\$us_os_win,
	"f|five_utr_gff:s" 	=>\$five_utr_gff,
	"g|genome_file:s" 	=>\$genome_file,
		
USAGE
}

if (!defined $put_ipt_gff){
	print_usage();
	exit;
}


## 1. Extend the 5'UTRs
system "bedtools slop -s -l 1 -r 0 -i $five_utr_gff -g $genome_file > temp.$$.1.gff";

## 2. Look for those putative IPTs that don't just overlap the 5'UTRs 
#(as they would in head to head genes where the 5'UTR of one gene has been called short)
system "bedtools intersect -s -v -a $put_ipt_gff -b temp.$$.1.gff > temp.$$.2.gff";


## 3. Find the gene id that the putative ipt is closest to (lost in previous script)
system "bedtools intersect -s -wa -wb -a temp.$$.2.gff -b $us_os_win > temp.$$.3.txt";

open (TMP, "<temp.$$.3.txt") or die "$!";
open (OUT, ">filtered.$put_ipt_gff") or die "$!";
while (<TMP>){
	chomp;
	my @g= split /\t/;
	print OUT "$g[0]\tPut_ipt_filters\tPutative_IPT\t$g[3]\t$g[4]\t$g[5]\t$g[6]\t$g[7]\t$g[17]\n";
}

exit;
