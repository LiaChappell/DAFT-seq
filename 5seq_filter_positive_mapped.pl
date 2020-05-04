#!/usr/bin/perl -w
use strict;

my $pos_sam =shift;
open (PAS, ">passed.$pos_sam") or die "$!";
open (MAYBE,">maybe.$pos_sam") or die "$!";
open (DEAD, ">failed.$pos_sam") or die "$!";
open (SAM, "<$pos_sam") or die "$!";
while (<SAM>){
	my @sam = split/\t/;
	my $cigar = $sam[5];
	my @soft = split (/S/, $cigar);
	if (($soft[0]>25) and ($soft[0] <35)){
		print PAS "$_";
	}
	elsif (($soft[0]>25) and ($soft[0] <40)) {
		print MAYBE "$_";
	}

	else {
		print DEAD "$_"; 
	}
}

