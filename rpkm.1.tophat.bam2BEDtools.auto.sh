#!/usr/local/bin/bash

if [ $# -ne 1 ]
then
    echo -e "\nUSAGE:  $0   <in_bam>\n"
    echo -e "DESCRIPTION: This script prepares a BAM file for subsequent intersection with a GFF file using BEDtools via Jason's other script tophat2.bam2RPKM_BEDtools.sh\n"
    echo -e "It does the following:"
    echo -e "* filter for minimum mapping quality of 30"
    echo -e "* adjust the CIGAR string of split reads (split reads spanning an intron), e.g.:"
    echo -e "  -  73M850N27M  is changed to  73M27S"
    echo -e "  - 40M2500N60M  is changed to  40S60M, plus the left-most mapping position POS of the read (4th column in the BAM file) is adjusted"
    echo -e "[ example: POS=1001 and CIGAR=40M2500N60M  are changed to  POS=3541 and CIGAR=40S60M ]"
    echo -e "\nNotes:"
    echo -e "* the length of the sequence reads is determined automatically"
    echo -e "* output is to file 'in_bam.forBedtoolsOnly.UNsorted.raw_mapQ30.bam'"
#    echo -e "* for a 2GB BAM file this script may take ~20min and ~1GB memory: with SORTing"
    echo -e "* for a 2GB BAM file this script may take ~10min and <100MB memory"
    echo -e "* the resulting BAM file will be UNsorted (some match starting positions are changed!)"
    echo -e "* the resulting BAM file will contain inconsistencies because some match starting positions are changed but neither corresponding read mates nor template lengths are changed! But this should not impact the GFF file intersection with BEDtools.\n"
    exit
fi

bam_in=$1
bam_out=$bam_in.forBEDToolsOnly.UNsorted.raw_mapQ30.bam

samtools view -h -q 30 $bam_in | perl -nle '@ar=split(/\t/); if ($ar[5] =~ /(\d+)M(\d+)N(\d+)M/){ $n1=$1;$n2=$2;$n3=$3; $readlen=length($ar[9]); if ($n1>=$n3) {$n3=($readlen-$n1);$ar[5]=$n1."M".($n3)."S"} else { $n1=($readlen-$n3);$ar[5]=($n1)."S".$n3."M";$ar[3]=$ar[3]+$n1+$n2;  } print join("\t",@ar)} else {print }'  | samtools view -Sb - > $bam_out

# Jason's original: samtools view -h -q 30 $bam_in | perl -nle '@ar=split(/\t/); if ($ar[5] =~ /(\d+)M(\d+)N(\d+)M/){ $n1=$1; $n3=$3; if ($n1>$n3) {$n3=(100-$n1);$ar[5]=$n1."M".($n3)."S"} else { $n1=(100-$n3);$ar[5]=($n1)."S".$n3."M"  } print join("\t",@ar)} else {print }'  | samtools view -Sb - > $bam_in.parsed.forBedtools.raw_mapQ30.bam
