#!/usr/bin/perl -w

# Reads in a VCF file and replaces missing genotypes with most frequent allele

use strict;
use warnings;
use List::Util qw(first min max);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
my $vcf_chr =$ARGV[1];
my $win_size =$ARGV[2];
my $win_step =$ARGV[3];



# ===========================================================================================
# Read in VCF file (requires AC and AN fields), record allele counts and calculate Hp in bins
# ===========================================================================================

# Open a file to output bin results to
my $file_bins ="$ARGV[0].dp";
open(my $fh_bins, '>', $file_bins) or die "Could not open file '$file_bins' $!";
#print $fh_bins "chr \t bin_start \t bin_stop \t sum_DP \t avg_DP \n";

# Initialize variables
my @win_pos;
my @win_bDP;
my $sum_DP =0;
my $avg_DP =0;
my $bin_start = 1;
my $bin_stop = $win_size;
my $tmp_start = 0;
my $tmp_stop = 0;
my $tmp_snps = 0;

# Open VCF file to read allele counts from
chomp($file);
open (FILE, $file); 

# Work through VCF file line by line
while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
		# Do nothing - ignore lines beginning with #
	} else {
		# Grab DP
		$_ = $data;
		my ($bDP) = /DP\=(.+?)[^0-9]/; 
                my @vcf = split(/\t/, $data);

		# Check values fall within bin before storing in array
		if ( $vcf[1] >= $bin_start && $vcf[1] <= $bin_stop )
		{
	                # Push data to array
	                push(@win_pos, $vcf[1]);
	                push(@win_bDP, $bDP);

		} else {

			# Reset allele counts and Hp
                 	$sum_DP = 0;
                 	$avg_DP = 0;

			# Calculate sum coverage for window
			for my $i ( reverse 0 .. $#win_pos) { $sum_DP = $sum_DP + $win_bDP[$i]; }
#			print "sum_DP $sum_DP \t win $#win_pos \n";
			if($sum_DP > 0) { $avg_DP = $sum_DP / (1+$#win_pos) } else { $avg_DP = 0 };

			# Move window forward 
			do {
				
				# Clear values if bin has no coverage, otherwise populate
				$tmp_start = 0;
				$tmp_stop = 0;
				if($win_pos[0]) { $tmp_start = $win_pos[0] };
				if($win_pos[-1]) { $tmp_stop = $win_pos[-1] };

				# Write vales to file
	                        print $fh_bins "$vcf_chr \t $bin_start \t $bin_stop \t $sum_DP \t $avg_DP \n";
	
				# Take a step along chromosome
				$bin_start += $win_step;
				$bin_stop += $win_step;
	
				# Remove elements from array that are out of bin dimensions
				for my $i (reverse 0 .. $#win_pos){
					if ( $win_pos[$i] < $bin_start || $win_pos[$i] > $bin_stop) {
	                                        splice @win_pos, $i, 1;
	                                        splice @win_bDP, $i, 1;
					}
				}
	
			} until ( $vcf[1] >= $bin_start && $vcf[1] <=  $bin_stop);
	
	                # Push new data to array
	                push(@win_pos, $vcf[1]);
	                push(@win_bDP, $bDP);

		}

	}
}


# Output last window, having reached end of file
for my $i ( reverse 0 .. $#win_pos) { $sum_DP = $sum_DP + $win_bDP[$i]; }
$avg_DP = $sum_DP / $#win_pos;

$tmp_start = 0;
$tmp_stop = 0;
if($win_pos[0]) { $tmp_start = $win_pos[0] };
if($win_pos[-1]) { $tmp_stop = $win_pos[-1] };
print $fh_bins "$vcf_chr \t $bin_start \t $bin_stop \t $sum_DP \t $avg_DP \n";



close FILE;
close $fh_bins;



exit;









