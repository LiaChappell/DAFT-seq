#! /usr/bin/perl -w
#
# File: spliceSites.bam2embl.pl
# Time-stamp: <27-Mar-2014 16:49:44 tdo>
# $Id: $
#
# Copyright (C) 2011 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#

my $limitAmount=0;

use strict;
use Data::Dumper;


### statics
my $minOverlap=shift;

###

my $emblDIR=shift;
my $bamfile=shift;

if (!defined($bamfile)) {
  print "Script to count confirm and potential new splice sites. This script looks for spliced reads in bam a file (N in cigar line), and reports their position. It does confirm known splice sites, and indicate possible alternative splicing.\n";
  print "\nusage: <anchor in bp - 1/3 of reads length> <directory with embl files> <any amount of bam files>  > Result.file.txt\n";
  print "\nLia advice: use 35 bp. Then list the split BAM files in same order as ls command. E.g. sample1.minus.bam sample1.plus.bam sample2.minus.bam sample2.plus.bam\n";
  
  exit
}

my @joined;
my @bamNames;

my ($h_knownSS,$h_product,$h_range,$h_strand,
	$h_sequence,$h_spliceL,$h_spliceR)=getEMBL($emblDIR);

while (defined($bamfile)) {
  my ($n)=$bamfile =~ /(\S+)\.bam/;
  
  push @bamNames, $n;
  print "work on $bamfile\n";
  
  my $h_newSplice = getSSBam($bamfile);
  push @joined, $h_newSplice;
  $bamfile=shift;
}



my %allSS;
my %allStrand;
my %newSpliceR;
my %newSpliceL;


# get all the splice sites into one hash
foreach my $i (0..(scalar(@joined)-1)) {
	foreach my $chr (keys %{ $joined[$i]}){
		foreach my $pos (keys %{ $joined[$i]{$chr}}){
			$allSS{$chr}{$pos}+=$joined[$i]{$chr}{$pos};
	#		print "$pos .... ";
			my ($first,$second)=$pos=~/^(\d+)\.\.(\d+)/;
			$newSpliceL{$chr}{$first}+=$joined[$i]{$chr}{$pos};
			$newSpliceR{$chr}{$second}+=$joined[$i]{$chr}{$pos};
			
			$allStrand{$chr}{$pos}=getStrand($h_sequence,$chr,$pos);
			
		}
	}
}

# include the one of the reference
foreach my $chr (keys %$h_knownSS){
	foreach my $pos (keys %{ $$h_knownSS{$chr}}){
		$allSS{$chr}{$pos}++;		
	}
}

### now print stuff
print "Chromosome\tStrand\tPosition\tGene";
foreach (@bamNames) {
  print "\t$_";
  
}
print "\tSum of reads per site\n";

my %stats;

foreach my $chr (sort keys %allSS){
  foreach my $pos (sort { ($a =~ /(\d+)\./)[0] <=> ($b =~ /(\d+)\./)[0]   }keys %{ $allSS{$chr}}){
	my $suml=0;  # counts the splice reads per site
	my ($first,$second)=$pos=~/^(\d+)\.\.(\d+)/;

	my $inGene=0;
	my $GeneStrand=-9;
	my $spliceStrand=-9;
	if (defined($allStrand{$chr}{$pos})) {
	  $spliceStrand=$allStrand{$chr}{$pos}
	} 
	
	my $possibleName=0;
	my $knownSS=0;

	my $leftinGene=0; my $rightinGene=0;
	
	
	
	if (!defined($$h_knownSS{$chr}{$pos})) {
	  print "$chr\t$allStrand{$chr}{$pos}\t$pos\t";
	}
	else {
	  print "$chr\t$$h_strand{$$h_knownSS{$chr}{$pos}}\t$pos\t";
	  $knownSS=1;
	}
	
	if (defined($$h_knownSS{$chr}{$pos})){
	  $inGene=1;
	  $GeneStrand=$$h_strand{$$h_knownSS{$chr}{$pos}};
	  print "$$h_knownSS{$chr}{$pos}"	
	} elsif(defined($$h_range{$chr}{$first})){
	  $inGene=1;
	  $GeneStrand=$$h_strand{$$h_range{$chr}{$first}};
	  print "$$h_range{$chr}{$first}"
	}
	elsif(defined($$h_range{$chr}{$second})) {
	  $inGene=1;
	  $GeneStrand=$$h_strand{$$h_range{$chr}{$second}};
	  
	  print "$$h_range{$chr}{$second}"
	}
	else {
	  $GeneStrand=-18;
	  
	  print "0"
		
	}
	
	### get which site is in gene
	if(defined($$h_range{$chr}{$first})){
	  $leftinGene=1;
	}
	if(defined($$h_range{$chr}{$second})) {
	  $rightinGene=1;
	}
	

	if ($knownSS){
	  $stats{$chr}{$pos}.="Known splice site";
	}
	else {
	   #defined($allStrand{$chr}{$pos})) {
	  
	  ### for alternative splicing, 
	  if ($inGene ) {
		
		if ($GeneStrand == $spliceStrand ) {
		  
		  ## check for exon skipping
		  if (!defined($$h_knownSS{$chr}{$pos}) &&
			  defined($$h_spliceL{$chr}{$first}) &&
			  defined($$h_spliceR{$chr}{$second})
			 ) {
			$stats{$chr}{$pos}.="Exon Skipping"
		  }
		  ### alternative start stpo
		  elsif ((!$leftinGene && $rightinGene) ||
				 ($leftinGene && !$rightinGene)
				) {
			$stats{$chr}{$pos}.="Alernative Start / Stop"
		  }
		  ### Altnative 3' 5'
		  elsif (
				 (defined($$h_spliceL{$chr}{$first}) && !defined($$h_spliceR{$chr}{$second}))
				 ||
				 (!defined($$h_spliceL{$chr}{$first}) && defined($$h_spliceR{$chr}{$second}))
				) {
			$stats{$chr}{$pos}.="3' 5' extension"
		  }
		  else {
			$stats{$chr}{$pos}.="two new"
		  }
		}# else some strand
		elsif ($spliceStrand!= -9) {
		  $stats{$chr}{$pos}.="Other Strand"
		} else {
		   $stats{$chr}{$pos}.="NonCanonical splice site"
		}
	  } # is part of a gene
	  else {
		$stats{$chr}{$pos}.="Not in Gene"
	  }
	}
	
	
	foreach my $i (0..(scalar(@joined)-1)) {
	  if (defined($joined[$i]{$chr}{$pos})){
		print "\t".$joined[$i]{$chr}{$pos};
		$suml+=$joined[$i]{$chr}{$pos}
	  } else {
		print "\t0"	
	  }
	}
	print "\t$suml\t";
	if (defined($stats{$chr}{$pos})) {
	  print $stats{$chr}{$pos}
	}
	if (defined($$h_knownSS{$chr}{$pos})) {
#if (defined($$h_product{$$h_knownSS{$chr}{$pos}})){	
	print "\t$$h_product{$$h_knownSS{$chr}{$pos}}";
 }
else { print "\t"};
	
	print "\n";

	### check for errors
	if ($knownSS && 0){
	  if ($suml==0 &&
		  ( ### one of the splice donor / accept set, report
		  defined($newSpliceL{$chr}{$first}) ||
		   defined($newSpliceR{$chr}{$second})
		   )
		 ) {
		if (defined($newSpliceL{$chr}{$first}) &&$newSpliceL{$chr}{$first}>5 ) {
		  print ">>> $chr $pos wrong splice site: check it... could be $first.. ??? $newSpliceL{$chr}{$first} \n"
		}
		if (defined($newSpliceR{$chr}{$second}) &&$newSpliceR{$chr}{$second}>5 ) {
		  print ">>> $chr $pos wrong splice site: check it... could be ?? .. $second ???$newSpliceR{$chr}{$second} \n"
		}
	  }
	}

	
  }
}


sub getStrand{
  my $h_seq = shift;
  my $chr   = shift;
  my $pos   = shift;

  my ($start,$stop) = $pos =~ /(\d+)\.\.(\d+)/;

 
  if (    uc(substr($$h_seq{$chr},$start,2))    eq "GT" &&
		  uc(substr($$h_seq{$chr},($stop-3),2)) eq "AG"
	 ) {
#	print "1 \n";
	
	return 1	
  }
  elsif ( uc(substr($$h_seq{$chr},$start,2)) eq "CT" &&
		  uc(substr($$h_seq{$chr},($stop-3),2)) eq "AC" ) {
#	print "0 \n";
	return 0
  }
  else {
#	print "-9\n";
	
	return -9
  }

}


sub getSSBam{
  my $file = shift;


  open F, "samtools view $file  | awk '\$5 > 5 && \$6 ~ \"N\"' |  cut -f 3,4,6 |" or die "Couldn't find bam file $file:$!\n";

  my %h;
  
  while (<F>) {
	my ($chr,$pos,$cigar) = split(/\t/);
	$chr =~ s/;//g;
	
	
	if ( $cigar =~ /(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M/) {
	#  print ;
	  $cigar =~ /(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M/;
	#  print "$pos $1 $2 $3 $4 $5\n"; 
	  my $start=$1+($pos-1);
	  my $stop=$2+(1+$start);
	#  print "$pos $1 $2 $3 $4 $5 -- $start -- $stop\n"; 
	  if ($1>$minOverlap || $3 >$minOverlap) {

		#print "set: $start .. $stop \n";
		$h{$chr}{"$start..$stop"}++;
	  }
	  $cigar =~ /(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M/;
	  if ($5 > $minOverlap || $3 >$minOverlap ) {
		$start=$stop+$3-1;
		$stop=$4+(1+$start);
#		print "$pos $1 $2 $3 $4 $5 -- $start -- $stop\n";
#		print "set2: $start .. $stop \n";
		$h{$chr}{"$start..$stop"}++;
	  }
	}
	elsif ($cigar =~ /(\d+)M(\d+)N(\d+)M/) {
	  
	  my ($start,$stop,$nextExon)= $cigar =~ /(\d+)M(\d+)N(\d+)M/;
	  #	print "$start $stop \n";
	  
	  ### we are just intested in "good aligments"
	  if ($start>$minOverlap && $nextExon > $minOverlap) {
		$start+=($pos-1);
		$stop+=(1+$start);
		#	print "$start $stop \n";
		
		$h{$chr}{"$start..$stop"}++;
	  }
	}
	
  }

  foreach my $chr (keys %h) {
	foreach my $pos (keys %{ $h{$chr}}) {
	  if ($h{$chr}{$pos}<=$limitAmount) {
		delete($h{$chr}{$pos})
	  }
	}
  }
  return (\%h)
	
}


sub getEMBL{
  my $dir =shift;
  
  opendir D, $dir or die "Problemes in open $dir: $! \n";
  my %h;
  my %product;
  my %range;
  my %Seq;
  my %Strand;
  my %spliceL;
  my %spliceR;
  
  map{
	my $file=$_;
	my $chr;
	
	if (/embl$/ || /embl.gz$/) {
	  if (/embl$/) {
		open F, "$dir/$file" or die "Problem open $file \n";
	  }
	  else {
		open F, "gunzip -c $dir/$file | " or die "Problem open $file \n";
	  }
	  ($chr)=$file =~ /^(\S+)\.embl/;
	  print $chr."\n";
	  
	  my @F=<F>;
#	  ($chr) = $F[0] =~ /ID\s+(\S+)\s;/;
	  
	  $chr =~ s/Transfer1\.//g;
	  $chr =~ s/\.final//g;
	   print $chr."\n";
	  close(F);

	  my $id;
	  my $product;
	  
	#  foreach (reverse(@F)){
	  foreach (my $i=((scalar(@F)-1)); $i>=0;$i--) {
		$_=$F[$i];
		chomp;
		if (/FT   CDS\s+(\S+)$/) {
### ATtention, know alternative splicing is ignoreed
		  if (! ($id=~/\.2$/ || $id=~/\.3$/ || $id=~/\.4$/)){
		  my $pos=$1;
		  $Strand{$id}=1;
		  if (/complement/) {
			$Strand{$id}=0;
		  }	
my $j=$i;
		  while ($pos =~ /,$/){
	#		$_=<>;
	#		print "pos= $pos\n$_\n";	
		$j++;
			
			$_=$F[$j];
#print "i = $i\n$_\n"; 
chomp;
			/FT   \s+(\S+)$/;
			$pos.=$1;
		  }
  #                print "$chr \t$pos\n";		  
		  my @ar=split(/\.\./,$pos);
		  $ar[0]=~/(\d+)$/;
		  my ($first)=$1;
		  $ar[(scalar(@ar)-1)]=~/^(\d+)/;
		  my ($last)	=$1;
		  
		  ### set the splice sites
		  foreach my $i (1..(scalar(@ar)-2)) {
			$ar[$i]=~/(\d+),(\d+)/;
			$h{$chr}{"$1..$2"}=$id;
			$spliceL{$chr}{$1}=$id;
			$spliceR{$chr}{$2}=$id;
			
			
		  }
		  for ($first..$last){
		  	$range{$chr}{$_}=$id;	
		  }

		  $Strand{$id}=1;
		  if (/complement/) {
			$Strand{$id}=0;
			
		  }
		} # end if alternative spliced
}
		elsif (/^FT  \s+\/product=\"(.*)$/) {
		  $product{$id}=$1;
		}
		elsif (/^FT  \s+\/systematic_id=\"(.*)\"/) {
			
		  $id=$1;
			$product{$id}="NA";
		} elsif (/^FT  \s+\/locus_tag=\"(.*)\"/) {
	#		print $id;
			  $id=$1;
      $product{$id}="NA";
		}
		
	  }

	  my $inSeq=0;
	  
	  ## now get the sequence
	  foreach (@F) {
	  
		if (/^SQ\s\s+Sequence\s+\d+/){
		  $inSeq=1;

		}elsif ($inSeq) {
		  
		  s/\d+//g; # awaz the number
		  s/\s//g;
		  s/\///g;
		  $Seq{$chr}.= $_;
		  
		}
	  }
	}
	
  }readdir (D);
  
  return (\%h, \%product, \%range, \%Strand, \%Seq, \%spliceL,\%spliceR);
  
}

