#!/usr/bin/perl -w
use strict;

my $neg_sam =shift;
open (PAS, ">passed.$neg_sam") or die "$!";
open (MAYBE,">maybe.$neg_sam") or die "$!";
open (DEAD, ">failed.$neg_sam") or die "$!";
open (SAM, "<$neg_sam") or die "$!";
while (<SAM>){
	my @sam = split/\t/;
	my $cigar = $sam[5];
	if ($cigar =~ /S$/){
		my @soft = split (/[A-Z]/, $cigar);
		if (($soft[-1]>25) and ($soft[-1] <35)){
			print PAS "$_";
		}
		elsif (($soft[0]>25) and ($soft[0] <40)) {
				print MAYBE "$_";
		}
	}
	else {
		print DEAD "$_"; 
	}
}
