#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my $us_gff;
my $cov_plus;
my $cov_minus;
my $output_file;

GetOptions
(
	"g|us_gff:s"					=> \$us_gff,
	"p|cov_plus:s"					=> \$cov_plus,
	"m|cov_minus:s"					=> \$cov_minus,
	"o|output_file:s"				=> \$output_file,
);

if (!defined $us_gff ||!defined $cov_plus ||!defined $cov_minus ||!defined $output_file ){
	print_usage();
	exit;
}

## File prep to use tabix

system "bgzip $cov_plus";
system "bgzip $cov_minus";
system "tabix -s 1 -b 2 -e 2 $cov_plus.gz";
system "tabix -s 1 -b 2 -e 2 $cov_minus.gz";



## Looking up coverage for each window

open (GFF, "<$us_gff" ) or die "$!";
open (OUT, ">$output_file") or die "$!";

while (<GFF>){
	chomp;
	my($z_chr, $z_source, $z_type, $z_start, $z_end, $z_score, $z_strand, $z_phase, $z_att) = split/\t/; #GFF
	my $max_cov_value = 0;
	my $max_cov_posit = 0;
	
	if ($z_strand eq "+"){
		open (TABIX_P, "tabix $cov_plus.gz $z_chr:$z_start-$z_end |") or die "$!";
		while (<TABIX_P>){
			chomp;
			my ($t_chr, $t_pos, $t_cov) = split/\t/;
			if ($t_cov > $max_cov_value) {
					$max_cov_value = $t_cov;
					$max_cov_posit = $t_pos;
			}
		}	
		close TABIX_P or die "$!"; # die important!
	}
	
	if ($z_strand eq "-"){
		open (TABIX_M, "tabix $cov_minus.gz $z_chr:$z_start-$z_end |") or die "$!";
		while (<TABIX_M>){
			chomp;
			my ($t_chr, $t_pos, $t_cov) = split/\t/;
			if ($t_cov > $max_cov_value) {
					$max_cov_value = $t_cov;
					$max_cov_posit = $t_pos;
			}
			
		}	
		close TABIX_M or die "$!"; # die important!
	}
		
	print OUT "$z_chr\tMax_TSS_Coverage_Finder\t5seq\t$max_cov_posit\t$max_cov_posit\t.\t$z_strand\t\.\t$z_att;Zone_Start:$z_start;Zone_End:$z_end;Max_Coverage=$max_cov_value;\n";
}

sub print_usage
{
	print <<USAGE;

	"g|us_gff:s"					=> \$us_gff,
	"p|cov_plus:s"					=> \$cov_plus,
	"m|cov_min:s"					=> \$cov_min,
	"o|output_file:s"				=> \$output_file,

USAGE
}
