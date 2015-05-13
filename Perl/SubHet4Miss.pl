#!/usr/bin/perl -w

# Substitutes heterozygous genotypes with missing genotypes

use strict;
use warnings;
use List::Util qw(first);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
chomp($file);
open (FILE, $file);
  

while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
		print "$data\n";
	} else {
		# Split VCF entry by tabs
		my @vcf = split(/\t/, $data);
		# Identify position of annotation for INFO field
		no warnings 'uninitialized';
		my $index = first { ((substr $vcf[$_], 0, 3) eq "GT:") } 0 .. $#vcf;

		my @INFO = split(/:/, $vcf[$index]);

                # Store all genotype data for the SNP in an array
                my @GTS = @vcf[$index+1 .. $#vcf];

                # Store alleles in an array
                my @alleles = map(substr($_,0,3),@GTS);

		# Substitute heterozyoug with missing
		my $imputed = $data;
#		$imputed =~ s/0\/1/\.\/\./g;
		$imputed =~ s/0\/1[^a-zA-Z0-9][0-9]*[^a-zA-Z0-9][0-9]*[^a-zA-Z0-9][0-9]*[^a-zA-Z0-9][0-9]*[^a-zA-Z0-9][0-9]*/\.\/./g;

# Original entry
$" = "\t";
#print "$data\n";
print "$imputed\n";


	}


}


close FILE;

exit;









