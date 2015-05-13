#!/usr/bin/perl -w

# Reads in a VCF file and replaces missing genotypes with most frequent allele

use strict;
use warnings;
use List::Util qw(first min max);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
my $nchr =$ARGV[1];
my $win_size =$ARGV[2];
my $win_step =$ARGV[3];



# ===========================================================================================
# Read in VCF file (requires AC and AN fields), record MAF in bins
# ===========================================================================================

# Open a file to output allele counts to
my $file_counts ="$ARGV[0].counts";
open(my $fh_counts, '>', $file_counts) or die "Could not open file '$file_counts' $!";
print $fh_counts "chr \t pos \t n_ref \t n_alt \t MAF \n";

# Open a file to output bin results to
my $file_bins ="$ARGV[0].hp";
open(my $fh_bins, '>', $file_bins) or die "Could not open file '$file_bins' $!";
print $fh_bins "bin_chr \t bin_start \t bin_stop \t chr \t prem_snp \t dern_snp \t n_snps \t sum_Maj \t sum_Min \t Hp \t MAF\n";

# Initialize variables
my @win_chr;
my @win_pos;
my @win_ref;
my @win_alt;
my @win_MAF;
my $sum_Maj =0;
my $sum_Min =0;
my $sum_MAF =0;
my $Hp = 0;
my $MAF = 0;
my $bin_chr = 1;
my $bin_start = 1;
my $bin_stop = $win_size;
my $tmp_chr = 0;
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
		# Grab AC and AN values to calculate #REF an #ALT counts
		$_ = $data;
		my ($AA) = /AC\=(.+?)[^0-9]/; 
		my ($AN) = /AN\=(.+?)[^0-9]/;
		my $AR = ($AN-$AA);
		my $MAF = ($AA/$AN);
		if ($MAF > 0.5) { ($MAF = (1-$MAF)); }
                my @vcf = split(/\t/, $data);
		print $fh_counts "$vcf[0] \t $vcf[1] \t $AR \t $AA \t $MAF\n";

		# Check values fall within bin before storing in array
		if ( $vcf[0] == $bin_chr && $vcf[1] >= $bin_start && $vcf[1] <= $bin_stop )
		{
	                # Push data to array
        	        push(@win_chr, $vcf[0]);
	                push(@win_pos, $vcf[1]);
	                push(@win_ref, $AR);
	                push(@win_alt, $AA);
			push(@win_MAF, $MAF);

		} else {

			# Reset allele counts and Hp
                 	$sum_Maj = 0;
                 	$sum_Min = 0;
			$sum_MAF = 0;
			$Hp = 0;
			$MAF = 0;

			# If current bin has coverage then output results
			if($#win_chr > 0){

				# Calculate sum of major and minor alleles
				for my $i ( reverse 0 .. $#win_pos)
				{
					if($win_ref[$i] >= $win_alt[$i]) { 
						$sum_Maj = $sum_Maj + $win_ref[$i];
						$sum_Min = $sum_Min + $win_alt[$i];
					} else {
	                                	$sum_Maj = $sum_Maj + $win_alt[$i];
	                               		$sum_Min = $sum_Min + $win_ref[$i];
					}
					$sum_MAF = $sum_MAF + $win_MAF[$i];
				}
				# Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
				if ($sum_Maj >0) { $Hp = (2 * $sum_Maj * $sum_Min) / ($sum_Maj + $sum_Min)**2; }
				if ($sum_MAF >0) { $MAF = ($sum_MAF / ($#win_chr +1)); }
			}

			# Move window forward 
			do {
				
				# Clear values if bin has no coverage, otherwise populate
				$tmp_chr = 0;
				$tmp_start = 0;
				$tmp_stop = 0;
				$tmp_snps = 0;
				if($win_chr[0]) { $tmp_chr = $win_chr[0] };
				if($win_pos[0]) { $tmp_start = $win_pos[0] };
				if($win_pos[-1]) { $tmp_stop = $win_pos[-1] };
				if($#win_chr) { $tmp_snps = ($#win_chr +1) };
				if($tmp_chr ==0) { 
					$Hp = 0;
					$MAF = 0;
					$sum_Maj = 0;
					$sum_Min = 0;
					$sum_MAF = 0;
					$tmp_snps = 0;
				}

				# Write vales to file
	                        print $fh_bins "$bin_chr \t $bin_start \t $bin_stop \t $tmp_chr \t $tmp_start \t $tmp_stop \t $tmp_snps \t $sum_Maj \t $sum_Min \t $Hp \t $MAF \n";
	
				# Take a step along chromosome
				$bin_start += $win_step;
				$bin_stop += $win_step;
	
				# Remove elements from array that are out of bin dimensions
				for my $i (reverse 0 .. $#win_pos){
					if ( $win_pos[$i] < $bin_start || $win_pos[$i] > $bin_stop) {
	                                        splice @win_chr, $i, 1;
	                                        splice @win_pos, $i, 1;
	                                        splice @win_ref, $i, 1;
	                                        splice @win_alt, $i, 1;
						splice @win_MAF, $i, 1;
					}
				}
	
				# If moving to different chromosome then reset bin and array
				if(  $vcf[0] != $bin_chr && $vcf[0] <= $nchr) {
					$bin_chr ++;
					$bin_start = 1;
					$bin_stop = $win_size;
					@win_chr =();
					@win_pos =();
					@win_ref =();
					@win_alt =();
					@win_MAF =();
				}
	
			} until ( $vcf[0] == $bin_chr && $vcf[1] >= $bin_start && $vcf[1] <=  $bin_stop);
	
	                # Push new data to array
	                push(@win_chr, $vcf[0]);
	                push(@win_pos, $vcf[1]);
	                push(@win_ref, $AR);
	                push(@win_alt, $AA);
			push(@win_MAF, $MAF);

		}

	}
}


# Output last window, having reached end of file
if($#win_chr > 0){
	# Calculate sum of major anbd minor alleles
        for my $i ( reverse 0 .. $#win_pos)
                {
                if($win_ref[$i] >= $win_alt[$i]) {
                        $sum_Maj = $sum_Maj + $win_ref[$i];
                	$sum_Min = $sum_Min + $win_alt[$i];
                } else {
                        $sum_Maj = $sum_Maj + $win_alt[$i];
                	$sum_Min = $sum_Min + $win_ref[$i];
        	}
		$sum_MAF = $sum_MAF + $win_MAF[$i];
        }
        # Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
	if ($sum_Maj >0) { $Hp = (2 * $sum_Maj * $sum_Min) / ($sum_Maj + $sum_Min)**2; }
	if ($sum_MAF >0) { $MAF = ($sum_MAF / ($#win_chr +1)); }
}
$tmp_chr = 0;
$tmp_start = 0;
$tmp_stop = 0;
$tmp_snps = 0;
if($win_chr[0]) { $tmp_chr = $win_chr[0] };
if($win_pos[0]) { $tmp_start = $win_pos[0] };
if($win_pos[-1]) { $tmp_stop = $win_pos[-1] };
if($#win_chr) { $tmp_snps = ($#win_chr +1) };
if($tmp_chr ==0) {
        $Hp = 0;
        $sum_Maj = 0;
        $sum_Min = 0;
        $tmp_snps = 0;
}
print $fh_bins "$bin_chr \t $bin_start \t $bin_stop \t $tmp_chr \t $tmp_start \t $tmp_stop \t $tmp_snps \t $sum_Maj \t $sum_Min \t $Hp \t $MAF \n";




close FILE;
close $fh_counts;
close $fh_bins;



exit;









