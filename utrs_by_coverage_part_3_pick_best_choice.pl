#!/usr/bin/perl
use strict;
use Getopt::Long;

##UTR combiner script##
my $gene_list;
my $utr_pref;
my $outpref;

GetOptions
(
	"l|gene_list:s"             => \$gene_list,
	"u|utr_pref:s"              => \$utr_pref,
);

if (!defined $utr_pref){
        print_usage();
        exit;
}

## Make arrays of 5' and 3' UTR files (assumes in order of preference!)
my @list_5utr =();
system "ls -v kept_5utrs.*.gff > list.$utr_pref.5utrs.txt";
open (LIS5, "<list.$utr_pref.5utrs.txt") or die "$!";
while (<LIS5>){
	chomp;
	push (@list_5utr, "$_");
}
	
my @list_3utr =();
system "ls -v kept_3utrs.*.gff > list.$utr_pref.3utrs.txt";
open (LIS3, "<list.$utr_pref.3utrs.txt") or die "$!";
while (<LIS3>){
    chomp;
	push (@list_3utr, "$_");
}


##Set up and "feed" 5UTR hash, where parent gene is key
my %kept_5utr =();
my %which_5utr_file =();

foreach (@list_5utr) {
	open(UTR5, "<$_") or die "$!";
	my $file_name = $_;
	while (<UTR5>){
        chomp;
        my ($data, $gene) = split/Parent=/;
		if (exists $kept_5utr{$gene}){
			#do nothing, have "best" 5UTR already
		}
		else{
			$kept_5utr{$gene} = $data;
			$which_5utr_file{$gene}=$file_name;
        }
	}
}

open (KEPT5, ">final.5utrs.$utr_pref.gff") or die "$!";
foreach my $gene (sort keys %kept_5utr){
	print KEPT5 "$kept_5utr{$gene}Parent=$gene\t$which_5utr_file{$gene}\n";
}


##Set up and "feed" 3UTR hash, where parent gene is key

my %kept_3utr =();
my %which_3utr_file =();

foreach (@list_3utr) {
    open(UTR3, "<$_") or die "$!";
	my $file_name = $_;
    while (<UTR3>){
        chomp;
        my ($data, $gene) = split/Parent=/;
        if (exists $kept_3utr{$gene}){
            #do nothing, have "best" 3UTR already
        }
        else{
            $kept_3utr{$gene} = $data;
            $which_3utr_file{$gene}=$file_name;
    	}
	}
}

open (KEPT3, ">final.3utrs.$utr_pref.gff") or die "$!";
foreach my $gene (sort keys %kept_3utr){
    print KEPT3 "$kept_3utr{$gene}Parent=$gene\t$which_3utr_file{$gene}\n";
}


##Read in list of all protein-coding genes

open (NO5, ">without.5utrs.$utr_pref.txt") or die "$!";
open (NO3, ">without.3utrs.$utr_pref.txt") or die "$!";

open (GENELIST, "<$gene_list") or die "$!";
while (<GENELIST>){
	chomp;
	if (!exists $kept_5utr{$_}){
		print NO5 "$_\tNo UTR detected using blocks in coverage method (or UTR failed one of the filters)\n";		
	}
	if (!exists $kept_3utr{$_}){
        	print NO3 "$_\tNo UTR detected using blocks in coverage method (or UTR failed one of the filters)\n";               
    	}
}

                


sub print_usage
{
	 print <<USAGE;
                
	"l|gene_list:s"             => Optional - will print genes that don't have UTR (just a simple text file, one id per line)
	"u|utr_pref:s"              => Prefix of your final UTR files - currently required...

USAGE
}
