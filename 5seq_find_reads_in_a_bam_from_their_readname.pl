#!/usr/bin/perl -w
use strict;

my $target_reads=shift;
my $inbam=shift;


my %read_1_keep;
my %read_2_keep;


open (READS, "<$target_reads") or die "$!"; 
while(<READS>){
	chomp;
	if (/1$/){
		s/\/1$//;
		$read_1_keep{$_} =1;
	}
	if (/2$/){
		s/\/2$//;
		$read_2_keep{$_} =2;
	}
#	{
#		print "Error with this input: $_\n";
#	}
}


open (F, "samtools view $inbam |") or die "$!";
open (KEEP, ">kept.$inbam.sam") or die "$!";

#my $readname;
#my $flag;

while (<F>){
	my ($readname, $flag) = split(/\t/);
	if ( (defined $read_1_keep{$readname}) and (0x0040 & $flag)){
		print KEEP "$_";
	}
	if ((defined $read_2_keep{$readname}) and (0x0080 & $flag)){
		print KEEP "$_";
	}
}

close F or die "$!";

#system(qq {samtools view $inbam | awk '\$1=="$_" && and(\$2,0x0040)'});
#system(qq {samtools view $inbam | awk '\$1=="$_" && and(\$2,0x0080)'});	
