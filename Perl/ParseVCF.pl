#!/usr/bin/perl -w

# Reads in a VCF file and extracts only homozygous variants
# No longer necessary, BCFtools performs this task admirably with ^het

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
		my @INFO = split(/:/, $vcf[-1]);
		if ($INFO[0] eq "0/0" || $INFO[0] eq "1/1") {
			$data =~ s/ /\t/g;
			print "$data\n";
		}
	}


}


close FILE;

exit;








