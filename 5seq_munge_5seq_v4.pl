#!/usr/bin/perl -w
use strict;

unless (@ARGV >0) {
        &USAGE;
}

sub USAGE {
	die "\n".'Usage: munge_5seq_v1.pl <input_bam> ## <genes_gff> <bedtools_genome_file>'."\n\n";
}

my $inbam = $ARGV[0];
my $genes=$ARGV[1];
my $genome_file=$ARGV[2];

## Other parameters
my $merge_dist =500;
my $sense_us_dist=5000;
my $sense_ds_dist=0;
my $antisense_us_dist=0;
my $antisense_ds_dist=100;


## Keep log file
#open (LOG, ">log.5tag.$$.txt") or die "$!";

## Make a file that just can the 5' end of TSS reads

# Convert BAM to BED, so can manipulate reads
system "bedtools bamtobed -i $inbam > 5tag.$$.bed";

# Reduce 5seq reads to 1 nt
# Account for strand when deciding which way to go 1 nt!
open (BED, "<5tag.$$.bed") or die "$!";
open (MOD_BED, ">mod_5tag.$$.bed") or die "$!";
while(<BED>){
	chomp;
	my($chr, $start, $end, $name, $score, $strand) = split/\t/; #BED
	#plus stand
	if ($strand eq "+"){
		my $mod_end = ($start +1);
		print MOD_BED "$chr\t$start\t$mod_end\t$name\t$score\t$strand\n";
	}
	#minus strand
	if ($strand eq "-"){
		my $mod_start = ($end -1);
		print MOD_BED "$chr\t$mod_start\t$end\t$name\t$score\t$strand\n";
	}  
}


## Get blocks of coverage for modified reads

# Get coverage for each base (per strand)
system "bedtools genomecov -i mod_5tag.$$.bed -g $genome_file -d -strand + > coverage.bga.plus.mod_5tag.$$.bed";
system "bedtools genomecov -i mod_5tag.$$.bed -g $genome_file -d -strand - > coverage.bga.minus.mod_5tag.$$.bed";

exit;

# Look for depth gtet 5 with 1nt reads (per strand)
#open (COV_P_1, "<coverage.bga.plus.mod_5tag.$$.bed") or die "$!";
#open (COV_P_1RES, ">gtet5.coverage.bga.plus.mod_5tag.$$.bed") or die "$!";
#while (<COV_P_1>){
#        my @covp= split /\s+/;
#        if ($covp[3] > 4){
#                print COV_P_1RES "$_";
#        }
#}
#close COV_P_1;
#close COV_P_1RES;

#open (COV_M_1, "<coverage.bga.minus.mod_5tag.$$.bed") or die "$!";
#open (COV_M_1RES, ">gtet5.coverage.bga.minus.mod_5tag.$$.bed") or die "$!";
#while (<COV_M_1>){
#        my @covm= split /\s+/;
#        if ($covm[3] > 4 ){
#                print COV_M_1RES "$_";
#        }
#}
#close COV_M_1;
#close COV_M_1RES;


# Merge regions with coverage gtet 5 into blocks (per strand)
#system "bedtools merge -d 100 -i gtet5.coverage.bga.plus.mod_5tag.$$.bed > merged.coverage.bga.plus.mod_5tag.$$.bed";
#system "bedtools merge -d 100 -i gtet5.coverage.bga.minus.mod_5tag.$$.bed > merged.coverage.bga.minus.mod_5tag.$$.bed";

# Make GFF files from BED files (to view in Artemis)
#system "cp merged.coverage.bga.plus.mod_5tag.$$.bed merged.coverage.bga.plus.mod_5tag.$$.gff";
#system "perl -p -i -e 's/v3/v3\tchado\t5seq/g' merged.coverage.bga.plus.mod_5tag.$$.gff";
#system "perl -p -i -e 's/\n/\t.\t+\t.\tID=\n/g' merged.coverage.bga.plus.mod_5tag.$$.gff";

#system "cp merged.coverage.bga.minus.mod_5tag.$$.bed merged.coverage.bga.minus.mod_5tag.$$.gff";
#system "perl -p -i -e 's/v3/v3\tchado\t5seqt/g' merged.coverage.bga.minus.mod_5tag.$$.gff";
#system "perl -p -i -e 's/\n/\t.\t-\t.\tID=\n/g' merged.coverage.bga.minus.mod_5tag.$$.gff";



## Get blocks of coverage for normal reads

# Get coverage for each base (per strand)
system "bedtools genomecov -i 5tag.$$.bed -g $genome_file -bga -strand + > coverage.bga.plus.5tag.$$.bed";
system "bedtools genomecov -i 5tag.$$.bed -g $genome_file -bga -strand - > coverage.bga.minus.5tag.$$.bed";


# Look at depth gtet 5 with normal reads (per strand)
open (COV_P, "<coverage.bga.plus.5tag.$$.bed") or die "$!";
open (COV_P_GTET, ">gtet5.coverage.bga.plus.5tag.$$.bed") or die "$!";
while (<COV_P>){
	my @covp= split /\s+/;
	if ($covp[3] > 4){
		print COV_P_GTET "$_";
	}
}
close COV_P;
close COV_P_GTET;

open (COV_M, "<coverage.bga.minus.5tag.$$.bed") or die "$!"; 
open (COV_M_GTET, ">gtet5.coverage.bga.minus.5tag.$$.bed") or die "$!";
while (<COV_M>){
	my @covm= split /\s+/;
	if ($covm[3] > 4 ){
		print COV_M_GTET "$_";
	}
}
close COV_M;
close COV_M_GTET;

# Merge regions with coverage gtet 5 into blocks (per strand)
system "bedtools merge -d $merge_dist -i gtet5.coverage.bga.plus.5tag.$$.bed > merged.coverage.bga.plus.5tag.$$.bed";
system "bedtools merge -d $merge_dist -i gtet5.coverage.bga.minus.5tag.$$.bed > merged.coverage.bga.minus.5tag.$$.bed";

# Make GFF files from BED files (to view in Artemis)
system "cp merged.coverage.bga.plus.5tag.$$.bed merged.coverage.bga.plus.5tag.$$.gff";
system "perl -p -i -e 's/v3/v3\tchado\t5seq/g' merged.coverage.bga.plus.5tag.$$.gff";
system "perl -p -i -e 's/\n/\t.\t+\t.\tID=\n/g' merged.coverage.bga.plus.5tag.$$.gff";

print "Hello!\n";
system "cp merged.coverage.bga.minus.5tag.$$.bed merged.coverage.bga.minus.5tag.$$.gff";
system "perl -p -i -e 's/v3/v3\tchado\t5seq/g' merged.coverage.bga.minus.5tag.$$.gff";
system "perl -p -i -e 's/\n/\t.\t-\t.\tID=\n/g' merged.coverage.bga.minus.5tag.$$.gff";
print "Hello again!\n";

## Filtering- keep region upstream of ORF

#Recap from top of script!
#my $sense_us_dist=2000;
#my $sense_ds_dist=0;
#my $antisense_us_dist=0;
#my $antisense_ds_dist=100;

print "Hello 1!\n";
# Get one file to work with
system "cat merged.coverage.bga.plus.5tag.$$.gff merged.coverage.bga.minus.5tag.$$.gff > pre_clean_tss_blocks_$$.gff";

# Get rid of noise from within genes
system "bedtools subtract -s -a pre_clean_tss_blocks_$$.gff -b $genes > post_clean_1_tss_blocks_$$.gff";

# Get rid of noise at 3' end of annotated genes (from their UTRs)
system "bedtools window -a $genes -b post_clean_1_tss_blocks_$$.gff -l $antisense_us_dist -r $antisense_ds_dist -sw -sm > probably_gene_end_noise_tss_blocks_$$.gff";
system "perl -p -i -e 's/^.*BST\t//g' probably_gene_end_noise_tss_blocks_$$.gff"; #Get rid of gene entries from the twice width gff file
system "perl -p -i -e 's/^.*GMT\t//g' probably_gene_end_noise_tss_blocks_$$.gff"; #Again

system "bedtools subtract -s -a post_clean_1_tss_blocks_$$.gff -b probably_gene_end_noise_tss_blocks_$$.gff > post_clean_2_tss_blocks_$$.gff";

# Now just get the region (usually 2kb) upstream of genes on the sense strand
system "bedtools window -a $genes -b post_clean_2_tss_blocks_$$.gff -l $sense_us_dist -r $sense_ds_dist -sw -sm > post_clean_3_tss_blocks_$$.gff";
system "perl -p -i -e 's/^.*BST\t//g' post_clean_3_tss_blocks_$$.gff"; #Get rid of gene entries from the twice width gff file
system "perl -p -i -e 's/^.*GMT\t//g' post_clean_3_tss_blocks_$$.gff"; #Again

## Read 1nt  coverage into hashes (for both strands)
# Bed files, so end coordinate is the same as start and end in gff
my %plus_coverage_1nt;
my %minus_coverage_1nt;

open (PLUS_COV, "<coverage.bga.plus.mod_5tag.$$.bed") or die "$!";
while (<PLUS_COV>){
	chomp;
	#my ($p_chr, $p_bed_start, $p_end, $p_cov) = split/\t/; #BED
	#my $p_position = "$p_chr:$p_end";
	my ($p_chr, $p_bed_start, $p_cov) = split/\t/; #BED
	my $p_start= ($p_bed_start + 1);
	my $p_position = "$p_chr:$p_start";
	#print "$p_position\n";
	#print "$p_cov\n";
	$plus_coverage_1nt{$p_position}=$p_cov;
}

open (MINUS_COV, "<coverage.bga.minus.mod_5tag.$$.bed") or die "$!";
while (<MINUS_COV>){
	chomp;
        #my ($m_chr, $m_bed_start, $m_end, $m_cov) = split/\t/; #BED
        #my $m_position = "$m_chr:$m_end";
	my ($m_chr, $m_bed_start, $m_cov) = split/\t/; #BED
	my $m_start= ($m_bed_start + 1);
        my $m_position = "$m_chr:$m_start";
	$minus_coverage_1nt{$m_position}=$m_cov;
}



## Now get max coverage for each window
# Use 1 nt read data

open (TSS_ZONES, "<post_clean_3_tss_blocks_$$.gff") or die "$!";
open (MAX_COV_POS, ">post_clean_4_tss_blocks_$$.gff") or die "$!";
while (<TSS_ZONES>){
	chomp;
	my $max_cov_for_zone = 0; 
	my $max_cov_position = 0;
	my($z_chr, $z_source, $z_type, $z_start, $z_end, $z_score, $z_strand, $z_phase, $z_att) = split/\t/; #GFF
	if ($z_strand eq "+"){
		for (my $i= $z_start; $i <= $z_end; $i++){
	 		print "$z_chr:$i\n";
			my $base_cov = $plus_coverage_1nt{"$z_chr:$i"};
			if ($base_cov > $max_cov_for_zone){
				$max_cov_for_zone = $base_cov;
				$max_cov_position = $i;
			}
		}
	}
	if ($z_strand eq "-"){
		for (my $i= $z_end; $i >= $z_start; $i++){
                        my $base_cov = $minus_coverage_1nt{"$z_chr:$i"};
                        if ($base_cov > $max_cov_for_zone){
                                $max_cov_for_zone = $base_cov;
                                $max_cov_position = $i;
                        }
                }
	}
	print MAX_COV_POS "$z_chr\tMax_TSS_Coverage_Finder\t5seq\t$max_cov_position\t$max_cov_position\t.\t$z_strand\t\.\tMax_Coverage=$max_cov_for_zone;Parent_id=\n";
}

exit;
