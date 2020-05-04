#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw( min max );

## File inputs


my $outpref; # prefix of coverage files, e.g. 3D7_v3
my $all_feat;  
my $parent_genes_gff;
my $parent_genes_rpkms;

my $thr;
my $parent_expr_thr=10;
my $olos_thr=20;


GetOptions
(
	"p|outpref:s"			=> \$outpref,
	"f|all_feat:s"                  => \$all_feat,
 	"g|parent_genes_gff:s"          => \$parent_genes_gff,
        "r|parent_genes_rpkms:s"        => \$parent_genes_rpkms,
	
	"t|thr:s"			=> \$thr,
	"e|parent_expr_thr:s"		=> \$parent_expr_thr,
	"o|olos_thr:s"			=> \$olos_thr,
);

if (!defined $parent_genes_gff ||!defined $parent_genes_rpkms ){
	print_usage();
	exit;
}

## Troubleshooting
#open (TEST, ">test.txt") or die "$!";

## Yell which inputs and threshold are in use
open (LOG, ">$outpref.thr$thr.log.txt") or die "$!";
print LOG "Inputs to script:\noutpref=ARGV[0]; # prefix of coverage files, e.g. 3D7_v3\nthr= ARGV[1]; # read coverage threshold, e.g. 5,10,15 .. 100\ngenes=ARGV[2]; #just the genes you want UTRs for, probably protein-coding genes\nall_feat=ARGV[3]; #include other annotated elements, such as ncRNAs etc\n\n";
print LOG "Starting analysis for a threshold of $thr reads\n";

## Get coverage at different ranges of depth- running this script once for each threshold
open (COV_P, "<$outpref.coverage.all.plus.strand.bed") or die "$!";
open (COV_P_GTET, ">$outpref.tmp.2.coverage.gtet.$thr.plus.strand.bed") or die "$!";
while (<COV_P>){
	my @covp= split /\s+/;
	if ($covp[3] > $thr){
		print COV_P_GTET "$_";
	}
}
close COV_P;
close COV_P_GTET;

open (COV_M, "<$outpref.coverage.all.minus.strand.bed") or die "$!"; 
open (COV_M_GTET, ">$outpref.tmp.2.coverage.gtet.$thr.minus.strand.bed") or die "$!";
while (<COV_M>){
	my @covm= split /\s+/;
	if ($covm[3] > $thr){
		print COV_M_GTET "$_";
	}
}
close COV_M;
close COV_M_GTET;

## Merge entries to give blocks of continuous coverage - "pseudo-molecules"
# Only positions which are next to each other- no gaps
# Merged bedfile has three columns: chr, feature start, feature end
system "bedtools merge -i $outpref.tmp.2.coverage.gtet.$thr.plus.strand.bed > $outpref.tmp.3.coverage.gtet.$thr.merged.plus.strand.bed";
system "bedtools merge -i $outpref.tmp.2.coverage.gtet.$thr.minus.strand.bed > $outpref.tmp.3.coverage.gtet.$thr.merged.minus.strand.bed";
print LOG "Merged coverage to make pseudo-molecules at read depth threshold of $thr\n";


## Load the blocks of coverage into an array for the next step, saves loading each transcript file many times!
open (TRANS_P, "<$outpref.tmp.3.coverage.gtet.$thr.merged.plus.strand.bed") or die "$!";
my @trans_p;
while(<TRANS_P>){
	chomp;
	push (@trans_p, $_);
}
open (TRANS_M, "<$outpref.tmp.3.coverage.gtet.$thr.merged.minus.strand.bed") or die "$!";
my @trans_m; 
while(<TRANS_M>){
	chomp;
	push (@trans_m, $_);
}

## Load in genes gff - this should only be starts and ends of genes- grep out the lines you want before this, no exons etc
open (GENES,"<$parent_genes_gff") or die "$!";

## Out files for next step
open (UTRSP, ">$outpref.tmp.4.coverage.gtet.$thr.utrs.plus.unfiltered.gff") or die "$!";
open (UTRSM, ">$outpref.tmp.4.coverage.gtet.$thr.utrs.minus.unfiltered.gff") or die "$!";


## For each gene, look for transcript (psuedo-molecules) that overlap it. Make the bits that overlap the end of the gene the UTRs
## Buyer beware, I'm looking for overlaps between a bed file (zero based indexing) and a gff file (one based indexing). My maths fudge should deal with that!
while(<GENES>){
	chomp;
	my($g_chr, $g_source, $g_type, $g_start, $g_end, $g_score, $g_strand, $g_phase, $g_att) = split/\t/; #GFF
	my(@att)= split (/;/, $g_att);
	my $gene_id = $att[0];
	$gene_id =~ s/^[a-zA-Z]+=//; #Tidies to just id, at least for GeneDB gffs!
	## This line might misbehave!

	#Plus strand transcripts (pseudomolecules)
	if ($g_strand eq "+"){
		foreach (@trans_p){
			my ($t_chr, $t_start_bed, $t_end)= split/\t/;
			my $t_start = ($t_start_bed +1); #bodge bed zero start to gff one start format
			if (($t_chr eq $g_chr) and ($t_start < $g_start) and ( $g_end < $t_end)){
				my $utr5_left = $t_start;
				my $utr5_right = ($g_start -1);
				my $utr3_left= ($g_end +1); 
				my $utr3_right= $t_end;
				print UTRSP "$g_chr\tUTR_finder_v2\t5UTR\t$utr5_left\t$utr5_right\t.\t+\t\.\tParent=$gene_id\n";
				print UTRSP "$g_chr\tUTR_finder_v2\t3UTR\t$utr3_left\t$utr3_right\t.\t+\t\.\tParent=$gene_id\n";
			}		
		} 
	}
	#Minus strand transcripts (pseudomolecules)
	if ($g_strand eq "-"){
		foreach (@trans_m){
			my ($t_chr, $t_start_bed, $t_end)= split/\t/;
			my $t_start = ($t_start_bed +1); #bodge bed zero start to gff one start format
			if (($t_chr eq $g_chr) and ($t_start < $g_start) and ( $g_end < $t_end)){
				my $utr5_left = ($g_end +1);
				my $utr5_right = $t_end;
				my $utr3_left= $t_start; 
				my $utr3_right= ($g_start -1);
				print UTRSM "$g_chr\tUTR_finder_v2\t5UTR\t$utr5_left\t$utr5_right\t.\t-\t\.\tParent=$gene_id\n";
				print UTRSM "$g_chr\tUTR_finder_v2\t3UTR\t$utr3_left\t$utr3_right\t.\t-\t\.\tParent=$gene_id\n";
			}
		} 
	}
}


## From the predicted UTRs, keep only those that don't overlap with a second gene (which implies that the psuedo-molecule at this depth threshold overlaps two annotated genes)
system "bedtools intersect -v -s -a $outpref.tmp.4.coverage.gtet.$thr.utrs.plus.unfiltered.gff -b $parent_genes_gff > $outpref.tmp.5.filtered.utrs.coverage.gtet.$thr.plus.strand.gff";
system "bedtools intersect -v -s -a $outpref.tmp.4.coverage.gtet.$thr.utrs.minus.unfiltered.gff -b $parent_genes_gff > $outpref.tmp.5.filtered.utrs.coverage.gtet.$thr.minus.strand.gff";
print LOG "Filtered out daft UTRs that overlap neighbouring genes on the SAME strand\n";

## Combine plus and minus strand UTRs for each depth threshold (keep 5UTR and 3UTR in same file for now)
system "cat $outpref.tmp.5.filtered.utrs.coverage.gtet.$thr.plus.strand.gff $outpref.tmp.5.filtered.utrs.coverage.gtet.$thr.minus.strand.gff > $outpref.tmp.6.filtered.unsorted.utrs.coverage.gtet.$thr.gff";
system "sort -k 1,1 -k 4,4n $outpref.tmp.6.filtered.unsorted.utrs.coverage.gtet.$thr.gff > $outpref.filtered.all.utrs.coverage.gtet.$thr.gff";
print LOG "Merged strands to give one gff file of UTRs per threshold";

## Get separate files for 5'UTRs and 3'UTRs
system "grep 5UTR $outpref.filtered.all.utrs.coverage.gtet.$thr.gff > $outpref.filtered.5utrs.coverage.gtet.$thr.gff";
system "grep 3UTR $outpref.filtered.all.utrs.coverage.gtet.$thr.gff > $outpref.filtered.3utrs.coverage.gtet.$thr.gff";

##Tidy clutter
print LOG "Check the word count of all the files just made before they are deleted\n";
my $wc = `wc *gtet.$thr.*`;
print LOG "\n\n$wc\n";
system "rm $outpref.tmp.*.$thr.*";
print LOG "Removed temporary files\n";

## Done with this threshold
print LOG "Done with munging for read threshold $thr\n";

########Part 2#############

## Fill Parent gene hashes (gene models, IDC RPKM)

# Hash = Parent gene gff (GeneDB format)
# Key = gene_id, value = rest of gff
my %parent_gff =();
open(PGFF,"<$parent_genes_gff") or die "$!";
while (<PGFF>) {
	chomp;
	my ($gff_bit, $gene_id)= split/ID=/; 
 	$gene_id =~ s/;.+\n/\n/;
	#$gene_id =~ s/^[a-zA-Z]+=//;
	$parent_gff{$gene_id} = $gff_bit;
}
close PGFF;

# Hash = Parent gene RPKMs (e,g, gene_id, 7 tps)
# Key = gene_id, value = max rpkm from time course
my %parent_max_rpkm =();
open (PRPKMS, "<$parent_genes_rpkms") or die "$!";
while (<PRPKMS>) {
	#print TEST "$_";
	chomp;
	my @tmp_prpkms = split/\t/;
	#print TEST "@tmp_prpkms\n";
	my $gene_id = shift @tmp_prpkms;
	#print TEST "$gene_id\n";
	my $max_rpkm = max @tmp_prpkms;
	#print TEST "Max=$max_rpkm.END\n";
	$parent_max_rpkm{$gene_id}= $max_rpkm;
}
close PRPKMS;

### 5UTRs ###
## Make a hash of UTR gff, where key is the gene id 
my %utr_5_gff = ();
open (THR_5GFF, "<$outpref.filtered.5utrs.coverage.gtet.$thr.gff") or die "$!";
while (<THR_5GFF>){
	chomp;
	#print TEST "$_\n";
	my ($gff_bits, $gene_id) = split/Parent=/;
		#$gene_id = s/^[a-zA-Z]+=//;
		#print TEST "Gene:$gene_id\n";
		$utr_5_gff{$gene_id}= $gff_bits;
	}
close THR_5GFF;
	
## Find UTRs that overlap genes on opposite strand - these are the UTRs that are likely to be antisense noise
system "bedtools intersect -wa -a $outpref.filtered.5utrs.coverage.gtet.$thr.gff -b $parent_genes_gff > tmp.5.$$.gff";
open (OLOS5, "<tmp.5.$$.gff") or die "$!";
# Fill hash with parent ids of UTRs that are opposite an annotated gene
# As before, key is gene_id - of the *UTR* parent
my %utr_5_olos = ();
while (<OLOS5>){
	chomp;
	my ($gff_bits, $gene_id) = split/Parent=/;
	$utr_5_olos{$gene_id}= $gff_bits;
	}
close OLOS5;
	
## Keep all UTRs that don't overlap genes on opposite strand (these are already cleaned in step 2)
## For UTRs that overlap genes on opp strand, keep only those where the max RPKM >50
## Keep only UTRs where parent gene is expressed gtet 10 RPKM
open(KEPT_5UTRS, ">kept_5utrs.$thr.gff") or die "$!";
open(FAIL_5UTRS, ">fail_5utrs.$thr.gff") or die "$!";

foreach my $utr_gene_id (sort keys %utr_5_gff){
	if((!exists $utr_5_olos{$utr_gene_id}) and ($parent_max_rpkm{$utr_gene_id} >= $parent_expr_thr)){
		print KEPT_5UTRS "$utr_5_gff{$utr_gene_id}Parent=$utr_gene_id\n";
	}
	elsif ((exists $utr_5_olos{$utr_gene_id}) and ($thr > $olos_thr) and ($parent_max_rpkm{$utr_gene_id} >= $parent_expr_thr)) {
		print KEPT_5UTRS "$utr_5_gff{$utr_gene_id}Parent=$utr_gene_id\n";
	}
	else{
		print FAIL_5UTRS "$utr_5_gff{$utr_gene_id}Parent=$utr_gene_id\n";
	}
}

## Make a hash of UTR gff, where key is the gene id 
my %utr_3_gff = ();
open (THR_3GFF, "<$outpref.filtered.3utrs.coverage.gtet.$thr.gff") or die "$!";
while (<THR_3GFF>){
        chomp;
        #print TEST "$_\n";
        my ($gff_bits, $gene_id) = split/Parent=/;
                #$gene_id = s/^[a-zA-Z]+=//;
                #print TEST "Gene:$gene_id\n";
                $utr_3_gff{$gene_id}= $gff_bits;
        }
close THR_3GFF;

## Find UTRs that overlap genes on opposite strand - these are the UTRs that are likely to be antisense noise
system "bedtools intersect -wa -a $outpref.filtered.3utrs.coverage.gtet.$thr.gff -b $parent_genes_gff > tmp.3.$$.gff";
open (OLOS3, "<tmp.3.$$.gff") or die "$!";
# Fill hash with parent ids of UTRs that are opposite an annotated gene
# As before, key is gene_id - of the *UTR* parent
my %utr_3_olos = ();
while (<OLOS3>){
        chomp;
        my ($gff_bits, $gene_id) = split/Parent=/;
        $utr_3_olos{$gene_id}= $gff_bits;
        }
close OLOS3;

## Keep all UTRs that don't overlap genes on opposite strand (these are already cleaned in step 2)
## For UTRs that overlap genes on opp strand, keep only those where the max RPKM >50
## Keep only UTRs where parent gene is expressed gtet 10 RPKM
open(KEPT_3UTRS, ">kept_3utrs.$thr.gff") or die "$!";
open(FAIL_3UTRS, ">fail_3utrs.$thr.gff") or die "$!";

foreach my $utr_gene_id (sort keys %utr_3_gff){
        if((!exists $utr_3_olos{$utr_gene_id}) and ($parent_max_rpkm{$utr_gene_id} >= $parent_expr_thr)){
                print KEPT_3UTRS "$utr_3_gff{$utr_gene_id}Parent=$utr_gene_id\n";
        }
        elsif ((exists $utr_3_olos{$utr_gene_id}) and ($thr > $olos_thr) and ($parent_max_rpkm{$utr_gene_id} >= $parent_expr_thr)) {
                print KEPT_3UTRS "$utr_3_gff{$utr_gene_id}Parent=$utr_gene_id\n";
        }
        else{
                print FAIL_3UTRS "$utr_3_gff{$utr_gene_id}Parent=$utr_gene_id\n";
        }
}


sub print_usage
{
	print <<USAGE;
	
	Input files: 
        "p|outpref:s"                   => Prefix of output from previous script- e.g. 3D7 
        "f|all_feat:s"                  => GFF with all features (that UTRs shouldn't overlap with on same strand)
        "g|parent_genes_gff:s"          => Genes that you want UTRs for- external coordinates, one line per gene
        "r|parent_genes_rpkms:s"        => Format gene, tab, number, tab etc - to only get UTRs from expressed genes
        
	Parameters: 
        "t|thr:s"                       => Threshold for coverage blocks **required**
        "e|parent_expr_thr:s"           => Threshold for calling parent genes expressed - default minimum is 10 RPKM
        "o|olos_thr:s"                  => Threshold for coverage blocks for calling overlapping UTRs - default is 20

USAGE
}




