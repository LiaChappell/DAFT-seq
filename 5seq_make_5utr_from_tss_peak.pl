#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $peaks_gff;
my $genes_gff; #external CDS coordinates
my $output_gff;

GetOptions
(
	"p|peaks_gff:s"					=> \$peaks_gff,
	"g|genes_gff:s"					=> \$genes_gff,	
	"o|output_gff:s"				=> \$output_gff,
);

if (!defined $peaks_gff ||!defined $genes_gff ||!defined $output_gff ){
	print_usage();
	exit;
}

#Get gene starts
open (GENES, "<$genes_gff" ) or die "$!";
my %gene_start; 
while (<GENES>){
	chomp;
	my($z_chr, $z_source, $z_type, $z_start, $z_end, $z_score, $z_strand, $z_phase, $z_att) = split/\t/; #GFF
	my $gene_id = $z_att;
	$gene_id =~ s/;.+$//g;
	$gene_id =~ s/^.+=//g;
	#print "$gene_id\n";
	if ($z_strand eq "+"){
		$gene_start{$gene_id}= $z_start;
	}

	if ($z_strand eq "-"){
		$gene_start{$gene_id}= $z_end;	
	}

}

#Get peak positions and print out matching gene starts
open (PEAK, "<$peaks_gff" ) or die "$!";
open (OUT, ">$output_gff") or die "$!";

while (<PEAK>){
        chomp;
        my($z_chr, $z_source, $z_type, $z_start, $z_end, $z_score, $z_strand, $z_phase, $z_att) = split/\t/; #GFF
        my $gene_id = $z_att;
        $gene_id =~ s/;.+$//;
        #$gene_id =~ s/^.+=//;
	#print "$gene_id\n";

        if ($z_strand eq "+"){
	my $utr_end = (($gene_start{$gene_id})-1); 
      	print OUT "$z_chr\t5UTR_based_on_TSO_3D7_IDC\t5UTR\t$z_start\t$utr_end\t.\t$z_strand\t\.\tParent=$gene_id\n";
        }

        if ($z_strand eq "-"){
        my $utr_end= (($gene_start{$gene_id}) + 1);
        print OUT "$z_chr\t5UTR_based_on_TSO_3D7_IDC\t5UTR\t$utr_end\t$z_start\t.\t$z_strand\t\.\tParent=$gene_id\n";
        }

}

sub print_usage
{
	print <<USAGE;

        "p|peaks_gff:s"                                 => \$peaks_gff,
        "g|genes_gff:s"                                 => \$genes_gff,
        "o|output_file:s"                               => \$output_file,

USAGE
}
