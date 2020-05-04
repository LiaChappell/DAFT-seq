#!/usr/bin/perl -w
use strict;
use Getopt::Long;

### GFF inputs ###
my $pseudo_mol; #make these with another script
my $all_feat; #gff with all features
my $utrs_5; #coverage-based 
my $utrs_3; #coverage-based
my $win; # putative ift windows
my $out; 

GetOptions
(
        "p|pseudo_mol:s"    		=> \$pseudo_mol,
	"a|all_feat:s"    		=> \$all_feat,
        "f|utrs_5:s"    		=> \$utrs_5,
	"t|utrs_3:s"    		=> \$utrs_3,
	"w|windows_of_put_ipts:s"	=> \$win,
	"o|out:s"   	 		=> \$out,
);

if (!defined $pseudo_mol || !defined $all_feat)
{
         print_usage();
         exit;
}

sub print_usage
{
        print <<USAGE;
        
	
	"p|pseudo_mol:s"       		=> \$pseudo_mol,
        "a|all_feat:s"          	=> \$all_feat,
        "f|utrs_5:s"            	=> \$utrs_5,
        "t|utrs_3:s"            	=> \$utrs_3,
        "w|windows_of_put_ipts:s"       => \$win,
	"o|out:s"               	=> \$out,               
        

USAGE
}


### 1. Round up all the things that IFTs shouldn't overlap with
system "cat $all_feat $utrs_5 $utrs_3 | sort -k 1,1 -k 4,4n > temp.$$.1.gff";


### 2. Keep the bits of the pseudo molecules that don't overlap with the other features
system "bedtools intersect -v -s -a $pseudo_mol -b temp.$$.1.gff > temp.$$.2.gff";


### 3. Look for the bits of pseudo molecules that are where you expect IPTs to be
system "bedtools intersect -wa -s -a temp.$$.2.gff -b $win > temp.$$.3.gff";

### 4. Merge all the fragments into windows to get something that's possible to work with
system "bedtools merge -s -d 100 -i temp.$$.3.gff > $out.gff";

exit;

 
 



