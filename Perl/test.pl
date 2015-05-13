#!/usr/bin/perl -w



use strict;
use warnings;
use List::Util qw(first);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
chomp($file);
open (FILE, $file);
  

sub most_frequent{ local *_=*_; $_[$_{$_}] = $_ for map{$_{$_}++; $_} +@_; $_[-1]; }

#print "Interval\t$chr\t$left\t$right\n";

while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
		print "$data\n";
	} else {
		# Split VCF entry by tabs
		my @vcf = split(/\t/, $data);
		# Identify position of annotation for INFO field
		my $index = first { ((substr $vcf[$_], 0, 3) eq "GT:") } 0 .. $#vcf;

		my @INFO = split(/:/, $vcf[$index]);

                # Store all genotype data for the SNP in an array
                my @GTS = @vcf[$index+1 .. $#vcf];

                # Store alleles in an array
                my @alleles = map(substr($_,0,3),@GTS);

                # Identify most frequent allele (Major allele)
                my $REF = grep (/0\/0/, @alleles);
                my $ALT = grep (/1\/1/, @alleles);
               	my $IMP = ();
                if ($REF > $ALT) { $IMP = "0/0" };
		if ($REF <= $ALT) { $IMP = "1/1" };

		# Impute Major allele
		my $imputed = $data;
		$imputed =~ s/\.\/\./$IMP/g;

# Original entry
$" = "\t";
#print "$data\n";
print "$imputed\n\n";


	}


}


close FILE;

exit;









