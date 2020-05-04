#!/bin/bash
set -e

if [ $# -ne 2 ]
then
   echo -e "\nusage: $0   <in.bam>   <outprefix>\n"
   exit
fi

inbam=$1
outpref=$2

samtools view -H $inbam > $outpref.plus.strand.sam
samtools view -F 4 $inbam | perl -nle 'print if /^[HM].+XS\:A\:\+/' >> $outpref.plus.strand.sam
samtools view -Sb $outpref.plus.strand.sam > $outpref.plus.strand.bam
samtools index $outpref.plus.strand.bam

samtools view -H $inbam > $outpref.minus.strand.sam
samtools view -F 4 $inbam | perl -nle 'print if /^[HM].+XS\:A\:\-/' >> $outpref.minus.strand.sam
samtools view -Sb $outpref.minus.strand.sam > $outpref.minus.strand.bam
samtools index $outpref.minus.strand.bam
