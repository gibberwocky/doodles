#!/usr/bin/perl -w

# Reads in a VCF file containing 2 individuals and creates a third which is
# a result of recombining allele 1 of individual 1 and allele 2 of individual 2
# ignoring missing sites

use strict;
use warnings;
use List::Util qw(first);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
my $id =$ARGV[1];
my $LN =0;

chomp($file);
open (FILE, $file);
  
while (my $data = <FILE>) {
	chomp $data;
	$data =~ s/AN=4/AN=6/g;

	if (substr($data, 0, 1) eq "#") {
		print "\n$data";
	} else {
		if ($LN == 0) {
		  $LN ++;
		  print "\t$id\n"; 
		}

		# Split VCF entry by tabs
		my @vcf = split(/\t/, $data);

		# Identify position of annotation for INFO field
		my $index = first { ((substr $vcf[$_], 0, 2) eq "GT") } 0 .. $#vcf;

		my @INFO = split(/:/, $vcf[$index]);

                # Store all genotype data for the SNP in an array
                my @GTS = @vcf[$index+1 .. $#vcf];

                # Store alleles in an array
                my @alleles = map(substr($_,0,3),@GTS);

		# Recombine non-missing genotypes
		if ($alleles[0] ne $alleles[1]) {
		  if ($alleles[0] ne './.') {
		    if ($alleles[1] ne './.') {
		      my $a = substr($alleles[0], 0, 1);
		      my $b = substr($alleles[1], 2, 1);
		      print "$data\t$a|$b:0:0.5,0,0.5\n";
		    }
		  }
		} else { print "$data\t$GTS[0]\n"; }

	}

}

close FILE;

exit;









