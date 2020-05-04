#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $genes_gff ;
my $genome_file;
my $us_dist;

GetOptions
(
	"g|genes_gff:s"					=> \$genes_gff,
	"f|genome_file:s"				=> \$genome_file,
	"u|us_dist:s"					=> \$us_dist,
);

if (!defined $genes_gff ||!defined $genome_file ||!defined $us_dist){
	print_usage();
	exit;
}

system "bedtools slop -s -l 0 -r 100 -i $genes_gff -g $genome_file > slopped_10bp_ds.$genes_gff";

system "bedtools flank -l $us_dist -r 0 -s -i $genes_gff -g $genome_file > unfiltered.$us_dist.bp_us.$genes_gff";

system "bedtools subtract -s -a unfiltered.$us_dist.bp_us.$genes_gff -b slopped_10bp_ds.$genes_gff > first_filter.$us_dist.of.$genes_gff";


## Tidy up bit of upstream regions that fall the other side of genes

# munge in genes gff
my %genes;
open (GENES, "<$genes_gff") or die "$!";
while (<GENES>){
	chomp;
	my($g_chr, $g_source, $g_type, $g_start, $g_end, $g_score, $g_strand, $g_phase, $g_att) = split/\t/;
	if ($g_strand eq "+"){
		$genes{$g_att} = $g_start; # 5' end of gene
	}
	if ($g_strand eq "-"){
		$genes{$g_att} = $g_end; # 5' end of gene
	}
}

# open zone gff, check that it's next to the gene it belongs to
open (KEEP, ">second_filter.$us_dist.of.$genes_gff") or die "$!";
open (Z_GFF, "<first_filter.$us_dist.of.$genes_gff") or die "$!";

while (<Z_GFF>){
	chomp;
	my($z_chr, $z_source, $z_type, $z_start, $z_end, $z_score, $z_strand, $z_phase, $z_att) = split/\t/; #GFF
	
	if ($z_strand eq "+"){
		if ( ($z_end +1) eq ($genes{$z_att}) ){ # if zone is actually next to gene
			print KEEP "$_\n";
		}
	}
	
	if ($z_strand eq "-"){
		if ( ($z_start -1) eq ($genes{$z_att}) ){ # if zone is actually next to gene
			print KEEP "$_\n";
		}
	}
}

# tidy up final field of gff
open (GFF, "<second_filter.$us_dist.of.$genes_gff") or die "$!";
open (FINAL, ">final_filtered.$us_dist.of.$genes_gff") or die "$!";

while(<GFF>){
	chomp;
	my($g_chr, $g_source, $g_type, $g_start, $g_end, $g_score, $g_strand, $g_phase, $g_att) = split/\t/; #GFF
	my(@att)= split (/;/, $g_att);
	my $gene_id = $att[0];
	$gene_id =~ s/^[a-zA-Z]+=//;
	print FINAL "$g_chr\t$g_source\t$g_type\t$g_start\t$g_end\t$g_score\t$g_strand\t$g_phase\t$gene_id\n";
}

system "rm unfiltered.$us_dist.bp_us.$genes_gff";
system "rm slopped_10bp_ds.$genes_gff";
system "rm first_filter.$us_dist.of.$genes_gff";
system "rm second_filter.$us_dist.of.$genes_gff";

sub print_usage
{
	print <<USAGE;
	
	"g|genes_gff:s"					=> \$genes_gff,
	"f|genome_file:s"				=> \$genome_file,
	"u|us_dist:s"					=> \$us_dist,	


USAGE
}

