#!/usr/bin/perl -w
use strict;
use Getopt::Long;

## Munges output from TDO's splicing script
## Assumes two bams for each sample: minus, plus

my $tdo_in; #output from TDO's script, minus and plus for each timepoint
my $out_pref; #prefix for out files
my $utrs_gff; #both 5' and 3' UTRs
my $genes_gff; #external coordinates of CDS only
my $sum_thr =5;


GetOptions
(
	"i|input_from_tdo_script:s"			=> \$tdo_in,
	"u|utrs_gff:s"					=> \$utrs_gff,
	"g|genes_gff:s"					=> \$genes_gff,
	"o|out_pref:s"					=> \$out_pref,
);

if (!defined $tdo_in ||!defined $out_pref || !defined $utrs_gff || !defined $genes_gff){
	print "\n$0 usage: -i tdo_input -u utrs_gff -g genes_gff -o out_pref\n\n";
	exit;
}

## Assorted outfiles
open (KNOWN, ">$out_pref.known_splice_sites.gff") or die "$!";
open (EXT, ">$out_pref.extension_within_genes.gff") or die "$!";
open (ALT_SS, ">$out_pref.alternative_start_stop_sites.gff") or die "$!";
open (EXON_SKIP, ">$out_pref.exon_skipping.gff") or die "$!";
open (OPPOSITE, ">$out_pref.opposite_genes.gff") or die "$!";
open (PUT_UTR_NCRNA, ">$out_pref.putative_utrs_or_ncrna.gff") or die "$!";
open (PUT_EXITRON_CAN, ">$out_pref.putative_exitrons_canonical.gff") or die "$!";
open (PUT_EXITRON_NON, ">$out_pref.putative_exitrons_noncanonical.gff") or die "$!";

## Munge infile
open (IN, "<$tdo_in") or die "$!"; 
my $seen_header=0;
my $header_line; 
while (<IN>){
	chomp;	
	if ($_ !~ /Chromosome/ && $seen_header==0){
		next;
	}
	if ($_ =~ /Chromosome/){
		$header_line = $_;
		$seen_header=1;
	}
	if ($_ !~ /Chromosome/ && $seen_header==1){
		my $ann ="NA";
		my @data = split /\t/;
		my $chr = shift @data;
		my $str = shift @data;
		my $pos = shift @data;
		my ($start, $end) = split (/\.\./, $pos);
		# Make it so that the splice site points don't include 1nt of CDS at each end. 
		$start= ($start +1);
		$end= ($end -1);
		my $id = shift @data;
		#print "$data[-2]\n";
		if ($data[-2] =~ /Known/){
			$ann = pop @data;
		}
		my $type = pop @data;
		my $sum = pop @data;
		my $array_size = @data;
		my $sum_plus=0;
		my $sum_minus=0;
		for (my $i=0; $i < ($array_size+1); $i +=2){
			$sum_minus= $sum_minus + $data[$i];
		}
		for (my $j=1; $j < ($array_size+1); $j +=2){
			$sum_plus= $sum_plus + $data[$j];
		}
		my $sum_both= $sum_plus +$sum_minus;
                my $guess_str="NA";
                if ($sum_plus > $sum_minus){
                        $guess_str="+";
                }
                if ($sum_minus > $sum_plus){
                        $guess_str="-";
                }
                if ($sum_plus eq $sum_minus){
                        $guess_str="unknown";
                }
                my $fixed_str;
                if ($str == 0){
                        $fixed_str="-";
                }
                if ($str == 1){
                        $fixed_str="+";
                }
                if ($str == -9){
                        $fixed_str= $guess_str;
                }

                ## GFF fields $chr, $source, $type, $start, $end, $score, $strand, $phase, $att
                my $gffline = "$chr\tSpliceSiteDetector\tSpliceSite\t$start\t$end\t.\t$fixed_str\t.\tID=$id;SumOfReads=$sum;SumBoth=$sum_both;SumOfReadsPlus=$sum_plus;SumOfReadsMinus=$sum_minus;GuessStrand=$guess_str\n";
		## GFF fields $chr, $source, $type, $start, $end, $score, $strand, $phase, $att
		#my $gffline = "$chr\tSpliceSiteDetector\tSpliceSite\t$start\t$end\t.\t$str\t.\tID=$id;SumOfReads=$sum;SumOfReadsPlus=$sum_plus;SumOfReadsMinus=$sum_minus\n";
		#print "Pos=$pos;Ann=$ann;Type=$type\n";
		#print "ID=$id;SUM=$sum\n";
		if ($sum >= $sum_thr){
			if ($type =~ /Known/){
				print KNOWN "$gffline";
			}
			if ($type =~ /extension/){
				print EXT "$gffline";
			}
			if ($type =~ /Alernative/){ #alt start stop
				print ALT_SS "$gffline";
			}
			if ($type =~ /Exon Skipping/){
				print EXON_SKIP "$gffline";
			}
			if ($type =~ /Other Strand/){
				print OPPOSITE "$gffline";
			}
			if ($type =~ /Not in Gene/){
				print PUT_UTR_NCRNA "$gffline";
			}
			if ($type =~ /two new/){
				print PUT_EXITRON_CAN "$gffline";
			}
			if ($type =~ /NonCanonical/){
				print PUT_EXITRON_NON "$gffline";
			}
		}
		
	}
}

####Munging GFFs
# $utrs_gff is both 5' and 3' UTRs

## 1. KNOWN
# Boring!

## 2. EXT
# Not much to say!

## 3. ALT_SS
# Find those that overlap with UTRs
system ("bedtools intersect -s -a $out_pref.alternative_start_stop_sites.gff -b $utrs_gff >$out_pref.munged.alternative_start_stop_sites.overlap_utrs.gff");
# Find those that overlaps with genes
system ("bedtools intersect -s -a $out_pref.alternative_start_stop_sites.gff -b $genes_gff >$out_pref.munged.alternative_start_stop_sites.overlap_utrs.gff");

## 4. EXON_SKIP
# Not much to say!

## 5. OPPOSITE
# Sort into those that do/don't overlap with UTRs
system ("bedtools intersect -s -a $out_pref.opposite_genes.gff -b $utrs_gff >$out_pref.munged.opposite_genes_ol_utrs.gff");
system ("bedtools intersect -v -s -a $out_pref.opposite_genes.gff -b $utrs_gff >$out_pref.munged.putative_ncrna_opp_genes.gff");

## 6. PUT_UTR_NCRNA
# Sort into those that do/don't overlap with UTRs
system ("bedtools intersect -s -a $out_pref.putative_utrs_or_ncrna.gff -b $utrs_gff >$out_pref.munged.spliced_utrs.gff");
system ("bedtools intersect -v -s -a $out_pref.putative_utrs_or_ncrna.gff -b $utrs_gff >$out_pref.munged.putative_ncrna_not_opp_genes.gff");

## 7. PUT_EXITRON_CAN
# Sort into those that overlap one gene/ two genes
# -c is for counts (final column)
system ("bedtools intersect -s -a $out_pref.putative_exitrons_canonical.gff -b $genes_gff -c > temp.$$.$out_pref.can_exitrons_ol_w_genes.gff");
open (TEMP_C, "<temp.$$.$out_pref.can_exitrons_ol_w_genes.gff") or die "$!";
open (CAN_0, ">$out_pref.munged.putative_exitrons_canonical.not_in_genes.gff") or die "$!";
open (CAN_1, ">$out_pref.munged.putative_exitrons_canonical.ol_1_gene.gff") or die "$!";
open (CAN_2, ">$out_pref.munged.putative_exitrons_canonical.ol_2_genes.gff") or die "$!";
 
while (<TEMP_C>){
	chomp;
	my @can = split /\t/;
	my $count= pop @can;
	if ($count =~ /0/){
		print CAN_0 join("\t", @can), "\n";
	}
	if ($count =~ /1/){
		print CAN_1 join("\t", @can), "\n";
	}
	if ($count =~ /2/){
		print CAN_2 join("\t", @can), "\n";
	}
}

close TEMP_C;
close CAN_0;
close CAN_1;
close CAN_2;

## 8. PUT_EXITRON_NON
# Sort into those that overlap one gene/ two genes
# -c is for counts (final column)
system ("bedtools intersect -s -a $out_pref.putative_exitrons_noncanonical.gff -b $genes_gff -c > temp.$$.$out_pref.noncan_exitrons_ol_w_genes.gff");
open (TEMP_N, "<temp.$$.$out_pref.noncan_exitrons_ol_w_genes.gff") or die "$!";
open (NON_0, ">$out_pref.munged.putative_exitrons_noncanonical.not_in_genes.gff") or die "$!";
open (NON_1, ">$out_pref.munged.putative_exitrons_noncanonical.ol_1_gene.gff") or die "$!";
open (NON_2, ">$out_pref.munged.putative_exitrons_noncanonical.ol_2_genes.gff") or die "$!";
 
while (<TEMP_N>){
	chomp;
	my @non = split /\t/;
	my $count= pop @non;
	if ($count =~ /0/){
		print NON_0 join("\t", @non), "\n";
	}
	if ($count =~ /1/){
		print NON_1 join("\t", @non),"\n";
	}
	if ($count =~ /2/){
		print NON_2 join("\t", @non),"\n";
	}
}

exit;
