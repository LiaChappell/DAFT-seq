#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;

my $help;
my $strand_values;
my $outfile; 

GetOptions
(
	"i|infile:s"	=> \$strand_values,
	"o|outfile:s"	=> \$outfile
);

if (!defined $strand_values)
{
		print_usage();
		exit;
}

my $fh =();
if (!defined $outfile)
{
		$fh = *STDOUT;
}
else
{
		open($fh, ">$outfile") or die "$!";
}


print $fh "1_gene_id\t2_length_nt\t3_strand\t4_sense_read_counts\t5_antisense_read_counts\t6_sense_rpkm\t7_antisense_rpkm\t8_sense_to_antisense_rpkm_ratio\n";

open(STRVAL, "<$strand_values") or die "$!: $strand_values";
my $header_line = <STRVAL>;

while(<STRVAL>)
{
	chomp;
	my($gene_id, $strand, $length_nt, $ps_read_counts, $ms_read_counts, $ps_rpkm, $ms_rpkm ) = split /\t/;
		

	if($strand eq "+" && $ms_rpkm ne "0")
	{
	my($plus_strand_sas_ratio) = sprintf("%.2f", $ps_rpkm/$ms_rpkm);
	print $fh "$gene_id"."\t"."$length_nt"."\t"."$strand"."\t"."$ps_read_counts"."\t"."$ms_read_counts"."\t"."$ps_rpkm"."\t"."$ms_rpkm"."\t"."$plus_strand_sas_ratio"."\n";
	}

	if($strand eq "+" && $ms_rpkm eq "0")
        {
        print $fh "$gene_id"."\t"."$length_nt"."\t"."$strand"."\t"."$ps_read_counts"."\t"."$ms_read_counts"."\t"."$ps_rpkm"."\t"."$ms_rpkm"."\t".$ms_rpkm."\n";
        }

	if($strand eq "-" && $ps_rpkm ne "0")
        {
	my($minus_strand_sas_ratio) = sprintf("%.2f", $ms_rpkm/$ps_rpkm);
        print $fh "$gene_id"."\t"."$length_nt"."\t"."$strand"."\t"."$ms_read_counts"."\t"."$ps_read_counts"."\t"."$ms_rpkm"."\t"."$ps_rpkm"."\t".$minus_strand_sas_ratio."\n";
	}

	if($strand eq "-" && $ps_rpkm eq "0")
	{
	print $fh "$gene_id"."\t"."$length_nt"."\t"."$strand"."\t"."$ms_read_counts"."\t"."$ps_read_counts"."\t"."$ms_rpkm"."\t"."$ps_rpkm"."\t"."\t".$ps_rpkm."\n";
	}
}
close STRVAL;

close $fh;

sub print_usage
{
print <<USAGE;

usage: rpkm.3.munge.table.to.sense.antisense.values.pl -i <infile> -o <outfile>
		
USAGE
}
