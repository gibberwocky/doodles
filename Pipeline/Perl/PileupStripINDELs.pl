#!/usr/bin/perl -w

# Reads in a VCF file and extracts non-INDEL variants

use strict;
use warnings;

my $file =$ARGV[0];
chomp($file);
open (FILE, $file);
  
#print "Interval\t$chr\t$left\t$right\n";

while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
		print "$data\n";
	} else {
		my @vcf = split(/\t/, $data);
		my @INFO = split(/;/, $vcf[7]);
		if ($INFO[0] ne 'INDEL') {
			local $" = "\t";
			print "@vcf\n";
		}
	}


}


close FILE;

exit;








#perl /home/plxdw1/Dropbox/Perl/PileupStripINDELs.pl 5677 JFM10_TAGCTT_L001_pileup_variants_hom.vcf

