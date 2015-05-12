#!/usr/bin/perl -w

use strict;
use warnings;

my $skip =$ARGV[0];
my $file =$ARGV[1];
my $sample =$ARGV[2];
chomp($file);
open (FILE, $file);
  
#print "Interval\t$chr\t$left\t$right\n";

while (my $data = <FILE>) {
	chomp $data;
	if ($. < $skip) {
		print "$data\n";
	} 
	if ($. == $skip) {
		chomp $sample;
		print"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
	}
	if ($. > $skip) {
		$data =~ s/ /\t/g;
		print "$data\n";
	} 


}


close FILE;

exit;








#perl /home/plxdw1/Dropbox/Perl/InsertVCFcolnames.pl 5677 ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_variants_SNPs_tmp.vcf JFM10_TAGCTT_L001 > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_variants_SNPs_hom.vcf

