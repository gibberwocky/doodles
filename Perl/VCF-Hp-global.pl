#!/usr/bin/perl -w

use strict;
use warnings;
use List::Util qw(first min max);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];


# ===========================================================================================
# Read in VCF file (requires AC and AN fields), record allele counts and calculate Hp
# ===========================================================================================

# Hp sum variables
my $sum_Maj =0;
my $sum_Min =0;
my $Hp = 0;

# Open VCF file to read allele counts from
chomp($file);
open (FILE, $file); 

# Work through VCF file line by line
while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
		# Do nothing - ignore lines beginning with #
	} else {
		# Grab AC and AN values to calculate #REF an #ALT counts
		$_ = $data;
		my ($AA) = /AC\=(.+?)[^0-9]/; 
		my ($AN) = /AN\=(.+?)[^0-9]/;
		my $AR = ($AN-$AA);

		# Record major and minor allele sums
		if ($AR>=$AA) {
                        $sum_Maj = $sum_Maj + $AR;
			$sum_Min = $sum_Min + $AA;
                } else {
                        $sum_Maj = $sum_Maj + $AA;
                        $sum_Min = $sum_Min + $AR;
                }

	}
}

$Hp = (2 * $sum_Maj * $sum_Min) / ($sum_Maj + $sum_Min)**2;

print "Genome-wide Hp = $Hp\n";


close FILE;



exit;









