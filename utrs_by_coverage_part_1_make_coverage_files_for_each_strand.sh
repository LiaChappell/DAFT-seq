#!/bin/bash

set -e

if [ $# -ne 2 ]
then
	echo -e "\nusage: $0 <in.bam> <outpref>\n"
	echo -e "UTRs by coverage (bedtools) and Adam's UTR finder script"
	echo -e "Need enough quota for 2 sam files from input bam!\n"
    exit
fi

### Lia's UTRs by coverage "pipeline" ###
# v1 (Sept 2014)
# for as much coverage as you can throw at the problem!
# designed for DAFT-seq (even coverage) in Pf 3D7 genome (compact, high AT content, overlapping UTRs)


## Pre-steps
# Map fastqs with TopHat v2, adding XS tag
# Merge BAM with Samtools, to get large BAM file (inbam)


## To run this script
# Run it in it's own directory
# Bsub with X Gb memory

## Inputs
inbam=$1
outpref=$2

## Split BAM into strands, using XS tag
echo -e "Let's start munging!"

samtools view -H $inbam > $outpref.plus.strand.sam
samtools view -H $inbam > $outpref.minus.strand.sam
echo -e "Made sam header for each strand"

samtools view -F 4 $inbam | perl -nle 'print if /^[HM].+XS\:A\:\+/' >> $outpref.plus.strand.sam
samtools view -F 4 $inbam | perl -nle 'print if /^[HM].+XS\:A\:\-/' >> $outpref.minus.strand.sam
echo -e "Munged bam into the two strands, based on XS tag. Also checked that read names start with H or M"

ls -lh *.sam > $outpref.sam.file.sizes.txt #for troubleshooting- all 4 exist on this line
echo -e "Listed the sam files for trouble shooting"

samtools view -Sb $outpref.plus.strand.sam > $outpref.plus.strand.bam
samtools view -Sb $outpref.minus.strand.sam > $outpref.minus.strand.bam
echo -e "Made bam for each strand"

rm $outpref.plus.strand.sam 
rm $outpref.minus.strand.sam 
echo -e "Removed filtered sam files"


## Get coverage at each position, separately by strand
# Note: I'm deliberately keeping the bits that are spliced, keeping the external coordinates of UTRs. Looking for spliced UTRs done later.
echo -e "Starting to calculate coverage"
bedtools genomecov -bga -ibam $outpref.plus.strand.bam > $outpref.coverage.all.plus.strand.bed
echo -e "Got plus strand coverage at each position"
bedtools genomecov -bga -ibam $outpref.minus.strand.bam > $outpref.coverage.all.minus.strand.bed
echo -e "Got minus strand coverage at each position"

#rm $outpref.plus.strand.bam
#rm $outpref.minus.strand.bam
#echo -e "Removed bam files"
echo -e "Leaving the plus and minus bam files for other scripts..."

echo -e "Only the two bedtools genome coverage files left now, other (large!) files removed!"
echo -e "UTRs_by_covarage_part_1 is done!"
