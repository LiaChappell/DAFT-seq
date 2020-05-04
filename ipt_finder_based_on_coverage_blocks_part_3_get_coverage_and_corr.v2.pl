#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::Basic qw(:all);
use List::Util qw(min max);

my $put_ipt_gff; #made by another script
my $five_utr_gff; #coverage-based 5'UTRs
my $bam_pref; #pref of bam files

GetOptions
(
	"i|put_ipt:s" 		=>\$put_ipt_gff,
	"f|five_utr_gff:s" 	=>\$five_utr_gff,
	"p|bam_pref:s"		=>\$bam_pref,
);

sub print_usage
{
	print <<USAGE;
		
	"i|put_ipt:s" 		=>\$put_ipt_gff,
	"f|five_utr_gff:s" 	=>\$five_utr_gff,
	"p|bam_pref:s"		=>\$bam_pref,
		
USAGE
}

if (!defined $put_ipt_gff || !defined $five_utr_gff || !defined $bam_pref){
	print_usage();
	exit;
}


## 1. List the BAMs files in the time course (assuming will ls in time order)

# List files
system "ls $bam_pref*plus.strand.bam > list_of_plus_bam_files.$$.txt";
system "ls $bam_pref*minus.strand.bam > list_of_minus_bam_files.$$.txt";

# Array to list each of the plus BAM files
my @plus_bam_files;
open (BAMP, "<list_of_plus_bam_files.$$.txt");
while(<BAMP>){
        chomp;
        push @plus_bam_files, $_;
        }

# Array to list each of the minus BAM files
my @minus_bam_files;
open (BAMM, "<list_of_minus_bam_files.$$.txt");
while(<BAMM>){
        chomp;
        push @minus_bam_files, $_;
        }

## 2. Bedtools coverage of 5'UTRs through the time course

# Loop through each of the plus BAM files
foreach my $loop_var (@plus_bam_files){
	system "bedtools coverage -sorted -abam $loop_var  -b $five_utr_gff > temp.$$.cov.$five_utr_gff.$loop_var.txt";
	#output: ($chr, $source, $type, $start, $end, $score, $strand, $phase, $att, $abs_cov, $bases_cov, $win_length, $frac_cov) #GFF #COV
	system "cut -f10 temp.$$.cov.$five_utr_gff.$loop_var.txt > temp.$$.useful_column.cov.$five_utr_gff.$loop_var.txt";
}

# Loop through each of the minus BAM files
foreach my $loop_var (@minus_bam_files){
	system "bedtools coverage -sorted -abam $loop_var  -b $five_utr_gff > temp.$$.cov.$five_utr_gff.$loop_var.txt";
	system "cut -f10 temp.$$.cov.$five_utr_gff.$loop_var.txt > temp.$$.useful_column.cov.$five_utr_gff.$loop_var.txt";
}

# Munge into useful tables
system "cut -f7 temp.$$.cov.$five_utr_gff.$plus_bam_files[0].txt > temp.$$.5utr_gene_ids_strand.txt";
system "cut -f9 temp.$$.cov.$five_utr_gff.$plus_bam_files[0].txt > temp.$$.5utr_gene_ids_plus.txt";
system "cut -f9 temp.$$.cov.$five_utr_gff.$minus_bam_files[0].txt > temp.$$.5utr_gene_ids_minus.txt";
system "paste temp.$$.5utr_gene_ids_plus.txt temp.$$.useful_column.cov.$five_utr_gff.*plus*.txt > temp.$$.5utr.coverage.plus.txt";
system "paste temp.$$.5utr_gene_ids_minus.txt temp.$$.useful_column.cov.$five_utr_gff.*minus*.txt > temp.$$.5utr.coverage.minus.txt";
system "paste temp.$$.5utr_gene_ids_strand.txt temp.$$.5utr.coverage.plus.txt temp.$$.5utr.coverage.minus.txt >temp.$$.5utr.coverage.both.txt";

# Make hashes to keep values in, key=gene id
my %five_utr_plus_cov_thr_idc;
my %five_utr_minus_cov_thr_idc;
my %five_utr_sense_cov_thr_idc;


# Populate the hashes
#open (PLUS_5, "<temp.13210.5utr.coverage.plus.txt") or die "$!";
open (PLUS_5, "<temp.$$.5utr.coverage.plus.txt") or die "$!";
while (<PLUS_5>){
	chomp;
	my ($gene_id, @cov_thru_idc)= split/\t/;
	my $vector = vector @cov_thru_idc;
	$five_utr_plus_cov_thr_idc{$gene_id} = $vector;
}

#open (MINUS_5, "<temp.13210.5utr.coverage.minus.txt") or die "$!";
open (MINUS_5, "<temp.$$.5utr.coverage.minus.txt") or die "$!"; 
while (<MINUS_5>){
	chomp;
	my ($gene_id, @cov_thru_idc)= split/\t/;
	my $vector = vector @cov_thru_idc;
	$five_utr_minus_cov_thr_idc{$gene_id} = $vector;
}

open (UTRS, "<$five_utr_gff") or die "$!";
open (CHECK, ">check.$$.5utr.coverage.txt") or die "$!";
while (<UTRS>){
	chomp;
	my($chr, $source, $type, $start, $end, $score, $strand, $phase, $att) = split/\t/; #GFF
	if ($strand eq "+"){
		$five_utr_sense_cov_thr_idc{$att}= $five_utr_plus_cov_thr_idc{$att};
	}
	if ($strand eq "-"){
		$five_utr_sense_cov_thr_idc{$att}= $five_utr_minus_cov_thr_idc{$att};
	}
	print CHECK "$att\t$strand\t$five_utr_sense_cov_thr_idc{$att}\t$five_utr_plus_cov_thr_idc{$att}\t$five_utr_minus_cov_thr_idc{$att}\n";
}

# Now in a mungable form!!


## 3. Bedtools coverage of putative IPTs through the time course

# Loop through each of the plus BAM files
foreach my $loop_var (@plus_bam_files){
	system "bedtools coverage -sorted -abam $loop_var  -b $put_ipt_gff > temp.$$.cov.$put_ipt_gff.$loop_var.txt";
	system "cut -f10 temp.$$.cov.$put_ipt_gff.$loop_var.txt > temp.$$.useful_column.cov.$put_ipt_gff.$loop_var.txt";
}

# Loop through each of the minus BAM files
foreach my $loop_var (@minus_bam_files){
	system "bedtools coverage -sorted -abam $loop_var  -b $put_ipt_gff > temp.$$.cov.$put_ipt_gff.$loop_var.txt";
	system "cut -f10 temp.$$.cov.$put_ipt_gff.$loop_var.txt > temp.$$.useful_column.cov.$put_ipt_gff.$loop_var.txt";
}

# Munge into useful tables
system "cut -f7 temp.$$.cov.$put_ipt_gff.$plus_bam_files[0].txt > temp.$$.put_ipt_gene_ids_strand.txt";
system "cut -f9 temp.$$.cov.$put_ipt_gff.$plus_bam_files[0].txt > temp.$$.put_ipt_gene_ids_plus.txt";
system "cut -f9 temp.$$.cov.$put_ipt_gff.$minus_bam_files[0].txt > temp.$$.put_ipt_gene_ids_minus.txt";
system "paste temp.$$.put_ipt_gene_ids_plus.txt temp.$$.useful_column.cov.$put_ipt_gff.*plus*.txt > temp.$$.put_ipt.coverage.plus.txt";
system "paste temp.$$.put_ipt_gene_ids_plus.txt temp.$$.useful_column.cov.$put_ipt_gff.*minus*.txt > temp.$$.put_ipt.coverage.minus.txt";
system "paste temp.$$.put_ipt_gene_ids_strand.txt temp.$$.put_ipt.coverage.plus.txt temp.$$.put_ipt.coverage.minus.txt > temp.$$.put_ipt.coverage.both.txt";



# Now in a mungable form!!
open (OUT, ">temp.$$.out.txt") or die "$!";
open (KEEP, ">temp.$$.keep.txt") or die "$!";
print OUT "gene_id\tipt_min_cov\tipt_max_cov\tcorr_5utr_ipt_sense\tcorr_ipt_sense_ipt_asens\t5utr_sense_cov\tipt_sense_cov\tipt_asens_cov\n";
print KEEP "gene_id\tipt_min_cov\tipt_max_cov\tcorr_5utr_ipt_sense\tcorr_ipt_sense_ipt_asens\t5utr_sense_cov\tipt_sense_cov\tipt_asens_cov\n";
open (IPT_COV, "<temp.$$.put_ipt.coverage.both.txt") or die "$!";
while (<IPT_COV>){
        #strand\tgene_id\tplus_cov\tgene_id\tminus_cov
	chomp;
        my @cov_both = split /\t/;
        my $strand= shift @cov_both;
        my $half = ((scalar @cov_both)/2);
        my @cov_plus = @cov_both[0..($half-1)];
        #print OUT "PLUS:@cov_plus\n";
        my @cov_minus = @cov_both[$half..$#cov_both];
        #print OUT "MINUS:@cov_minus\n";
	my $gene_id = shift @cov_plus;
	shift @cov_minus; #throw away 2nd gene id
	my @cov_sense = ();
	my @cov_asens = ();
        if ($strand eq "+"){
 	        #print OUT "Strand=$strand; PLUS:@cov_plus\n";
		@cov_sense = @cov_plus;
		@cov_asens = @cov_minus
        }
        if ($strand eq "-"){
                #print OUT "Strand=$strand; MINUS:@cov_minus\n";
		@cov_sense = @cov_minus;
		@cov_asens = @cov_plus;
        }

	my $ipt_cov_min = min @cov_sense;
	my $ipt_cov_max = max @cov_sense;

	my $vector_ipt_sense = vector @cov_sense;
	my $vector_ipt_asens = vector @cov_asens;
	my $vector_5utr_sense = ($five_utr_sense_cov_thr_idc{$gene_id});

        my $corr_5utr_ipt_sense= correlation ($vector_5utr_sense, $vector_ipt_sense);
        my $corr_ipt_sense_ipt_asens = correlation ($vector_ipt_sense, $vector_ipt_asens);
        print OUT "$gene_id\t$ipt_cov_min\t$ipt_cov_max\t$corr_5utr_ipt_sense\t$corr_ipt_sense_ipt_asens\t$vector_5utr_sense\t$vector_ipt_sense\t$vector_ipt_asens\n";
	if ($ipt_cov_max>4 && $corr_5utr_ipt_sense>0.6 && $corr_ipt_sense_ipt_asens<0.5){
		print KEEP "$gene_id\t$ipt_cov_min\t$ipt_cov_max\t$corr_5utr_ipt_sense\t$corr_ipt_sense_ipt_asens\t$vector_5utr_sense\t$vector_ipt_sense\t$vector_ipt_asens\n";
	}
}

exit;
