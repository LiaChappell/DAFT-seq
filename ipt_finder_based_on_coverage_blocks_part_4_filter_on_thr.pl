#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $put_ipt_file;
my $min_cov;
my $max_cov;
my $five_utr_corr;
my $as_corr;
my $h2h_ids;

GetOptions
(
	"p|put_ipt_file:s"	=>\$put_ipt_file,
	"l|min_cov:s"		=>\$min_cov,
	"h|max_cov:s"		=>\$max_cov,
	"f|five_utr_corr:s"	=>\$five_utr_corr,
	"a|as_corr:s"		=>\$as_corr,
	"i|head_2_heads_ids:s"	=>\$h2h_ids,
);

if (!defined $put_ipt_file){
	print_usage();
	exit;
}

sub print_usage
{
	print <<USAGE;

        "p|put_ipt_file:s"      =>\$put_ipt_file,
        "l|min_cov:s"           =>\$min_cov,
        "h|max_cov:s"           =>\$max_cov,
        "f|five_utr_corr:s"     =>\$five_utr_corr,
        "a|as_corr:s"           =>\$as_corr,
	"i|head_2_heads_ids:s"  =>\$h2h_ids,
	
USAGE
}

my %h2h_hash;
open (H2H, "<$h2h_ids") or die "$!";
while (<H2H>){
	chomp;
	$h2h_hash{$_}=1;
}


open (IN, "<$put_ipt_file") or die "$!";
open (OUT, ">filtered.$put_ipt_file.min_$min_cov.max_$max_cov.five_utr_corr_$five_utr_corr.as_corr_$as_corr.txt") or die "$!";
while (<IN>){
	chomp;
	my @ipt = split /\t/;
	#GFF_1	GFF_2	GFF_3	GFF_4	GFF_5	GFF_6	GFF_7	GFF_8	GFF_9	gene_id	ipt_min_cov	ipt_max_cov	corr_5utr_ipt_sense	corr_ipt_sense_ipt_asen5utr_sense_cov	ipt_sense_cov	ipt_asens_cov
	my $gene_id =$ipt[9];
	$gene_id =~ s/Parent=//;
	$gene_id =~ s/;.+//;
	print "$gene_id\n";
	my $l_min_cov =$ipt[10];
	my $l_max_cov  =$ipt[11];
	my $l_five_utr_corr  =$ipt[12];
	my $l_as_corr  =$ipt[13];
	if ((!exists $h2h_hash{$gene_id}) && ($l_min_cov >= $min_cov) && ($l_max_cov >= $max_cov) && ($l_five_utr_corr >=$five_utr_corr) && ($as_corr >= $l_as_corr) && ($l_five_utr_corr >= $l_as_corr)){
		print OUT "$_\n";	
	}
}
